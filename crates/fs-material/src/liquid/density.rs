//! Local liquid thermodynamics from the original two-coordinate density law.
//!
//! The selected `rho(T, p)` source and the same liquid state's Cp determine
//! thermal expansion, isothermal/adiabatic compressibility, Cv and sound speed.
//! Derivatives are analytical tangents of the declared affine source cells;
//! no perturbed queries, fitted EOS, extrapolation or independent sound-speed
//! claim is introduced. The returned liquid retains its original property
//! receipts. The derived response remains separate from that source evidence.
//!
//! This is a local, single-phase, fixed-composition equilibrium response. It
//! does not establish phase stability, derivative uncertainty, or integrability
//! of independently sampled Cp and density into a global thermodynamic potential.

use core::fmt;
use std::sync::Arc;

use fs_matdb::{
    ClaimId, InterpolationPolicy, MaterialCard, PreparedSurfaceResponse, PropertyValue, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermochem::MassSpecificThermalQuantityV1;
use fs_thermochem::eos::density::{
    DensityPressureDerivativeV1, DensityResponseErrorV1, DensityTemperatureDerivativeV1,
    DensityThermodynamicResponseV1,
};

use super::ResolvedLiquidState;
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};

/// Original source, thermodynamic-stability or shared fluid-view refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum LiquidDensityError {
    /// The original card, context, density source or property refresh refused.
    Source(MaterialStatePointError),
    /// The source derivatives and Cp do not admit a finite stable local response.
    Response(DensityResponseErrorV1),
    /// A derived acoustic/transport coefficient is not representable.
    Fluid(MaterialError),
}

impl fmt::Display for LiquidDensityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Source(source) => write!(f, "liquid density source: {source}"),
            Self::Response(source) => source.fmt(f),
            Self::Fluid(source) => source.fmt(f),
        }
    }
}
impl std::error::Error for LiquidDensityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Source(source) => source,
            Self::Response(source) => source,
            Self::Fluid(source) => source,
        })
    }
}
impl From<MaterialStatePointError> for LiquidDensityError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::Source(source)
    }
}
impl From<DensityResponseErrorV1> for LiquidDensityError {
    fn from(source: DensityResponseErrorV1) -> Self {
        Self::Response(source)
    }
}
impl From<MaterialError> for LiquidDensityError {
    fn from(source: MaterialError) -> Self {
        Self::Fluid(source)
    }
}

#[derive(Clone, Debug)]
struct DensitySource {
    surface: PreparedSurfaceResponse,
    context: QueryPoint,
    property: String,
    claim: ClaimId,
    temperature_index: usize,
    validity: [Option<(f64, f64)>; 2],
}

/// The original liquid with a coherent local caloric/acoustic density response.
/// Immutable prepared source data is shared across all refreshed states.
#[derive(Clone, Debug)]
pub struct ResolvedLiquidDensityState {
    liquid: ResolvedLiquidState,
    source: Arc<DensitySource>,
    response: DensityThermodynamicResponseV1,
    thermoviscous: ThermoviscousFluidState,
    thermal: ThermalFluidTransportState,
}

impl ResolvedLiquidDensityState {
    /// Bind the liquid's already selected density claim as `rho(T, p)`.
    ///
    /// Both named axes require their exact source dimensions and an explicitly
    /// typed query: absolute-temperature/static-pressure semantics or matching
    /// dimension-only coordinates. Only original `LinearInside` Surface2D
    /// cells with a unique two-sided gradient qualify. Source edges, holes,
    /// unequal incident gradients and declared T/p validity endpoints refuse.
    /// A common gradient across an interior cell boundary remains admissible.
    ///
    /// An existing density-expansion path or independently selected beta refuses:
    /// this owner never replaces another thermal-expansion authority silently.
    /// Original source caps bound preparation and subsequent gradient work.
    pub fn try_from_liquid(
        card: &MaterialCard,
        liquid: ResolvedLiquidState,
        temperature_axis: &str,
        pressure_axis: &str,
    ) -> Result<Self, LiquidDensityError> {
        if liquid.material().card_identity() != card.content_hash() {
            return Err(invalid(
                "liquid density response requires the original immutable card",
            ));
        }
        if liquid.volumetric_expansion_per_k().is_some() {
            return Err(invalid(
                "liquid density response cannot replace an existing expansion source",
            ));
        }
        let context = liquid.material().query();
        validate_coordinates(&context, temperature_axis, pressure_axis)?;
        let property = liquid.requirements().properties()[0].name().to_owned();
        let selected = liquid
            .material()
            .property(&property)
            .expect("liquid retains its required density")
            .answer()
            .receipt
            .selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or_else(|| invalid("selected liquid density claim is absent"))?;
        let PropertyValue::Surface2D {
            axes, axis_dims, ..
        } = &claim.value
        else {
            return Err(invalid(
                "liquid density response needs a two-coordinate source surface",
            ));
        };
        let temperature_index = axes
            .iter()
            .position(|axis| axis == temperature_axis)
            .ok_or_else(|| invalid("liquid density source has no selected temperature axis"))?;
        if axes[1 - temperature_index] != pressure_axis
            || axis_dims[temperature_index] != Temperature::DIMS
            || axis_dims[1 - temperature_index] != Pressure::DIMS
            || claim.interpolation != InterpolationPolicy::LinearInside
        {
            return Err(invalid(
                "liquid density source needs original linear T/p cells with matching dimensions",
            ));
        }
        let surface = claim.value.prepare_surface_response().map_err(|source| {
            MaterialStatePointError::Query {
                property: property.clone(),
                source,
            }
        })?;
        let source = Arc::new(DensitySource {
            surface,
            context,
            property,
            claim: selected,
            temperature_index,
            validity: axes.each_ref().map(|axis| claim.validity.bound(axis)),
        });
        Self::sample(liquid, source)
    }

    /// Refresh every original liquid property and the response at the SAME T/p.
    /// The original card, claims, coordinate meanings and all other values
    /// (including composition and shear rate) remain fixed. A failure returns
    /// no partially refreshed state; source tables are shared, not recompiled.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, LiquidDensityError> {
        self.source.require_context(point)?;
        Self::sample(
            self.liquid.at_query_point(card, point)?,
            Arc::clone(&self.source),
        )
    }

    /// Original liquid transport state, without added or manufactured claims.
    pub const fn liquid(&self) -> &ResolvedLiquidState {
        &self.liquid
    }
    /// Complete original source bundle and usage receipts at this state point.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.liquid.material()
    }
    /// Derived local response; this does not alter the original material bundle.
    pub const fn response(&self) -> &DensityThermodynamicResponseV1 {
        &self.response
    }
    /// Explicit original source coordinate carrying absolute temperature.
    pub fn temperature_axis(&self) -> &str {
        &self.source.surface.axes()[self.source.temperature_index]
    }
    /// Explicit original source coordinate carrying thermodynamic pressure.
    pub fn pressure_axis(&self) -> &str {
        &self.source.surface.axes()[1 - self.source.temperature_index]
    }

    fn sample(
        liquid: ResolvedLiquidState,
        source: Arc<DensitySource>,
    ) -> Result<Self, LiquidDensityError> {
        let point = liquid.material().query();
        source.require_context(&point)?;
        let property = liquid
            .material()
            .property(&source.property)
            .ok_or_else(|| invalid("liquid density refresh lost its original source"))?;
        if property.answer().receipt.selected != source.claim {
            return Err(invalid("liquid density refresh changed its selected claim"));
        }
        let coordinates = source
            .surface
            .axes()
            .each_ref()
            .map(|axis| point.axes()[axis]);
        for (value, bounds) in coordinates.into_iter().zip(source.validity) {
            if bounds.is_some_and(|(lo, hi)| value <= lo || value >= hi) {
                return Err(invalid(
                    "liquid density derivative needs strict interior T/p validity",
                ));
            }
        }
        let gradient = source
            .surface
            .gradient_for_axes(coordinates, [true; 2])
            .map_err(|source_error| MaterialStatePointError::Query {
                property: source.property.clone(),
                source: source_error,
            })?;
        // The prepared cells belong to this exact selected claim. Keep its
        // canonical queried density rather than rebuilding it from a tangent.
        let response = DensityThermodynamicResponseV1::try_new(
            Temperature::new(coordinates[source.temperature_index]),
            Density::new(liquid.density_kg_m3()),
            MassSpecificThermalQuantityV1::new(liquid.specific_heat_j_kg_k()),
            DensityTemperatureDerivativeV1::new(gradient[source.temperature_index]),
            DensityPressureDerivativeV1::new(gradient[1 - source.temperature_index]),
        )?;
        let transport = liquid.transport_state()?;
        let acoustic =
            FluidAcousticState::try_new(liquid.density_kg_m3(), response.sound_speed_m_s())?;
        let thermoviscous = ThermoviscousFluidState::try_new(
            acoustic,
            transport,
            response.specific_heat_cv_j_kg_k(),
        )?;
        let thermal = ThermalFluidTransportState::try_new(
            transport,
            response.temperature_k(),
            Some(response.volumetric_expansion_per_k()),
        )?;
        Ok(Self {
            liquid,
            source,
            response,
            thermoviscous,
            thermal,
        })
    }
}

impl DensitySource {
    fn require_context(&self, point: &QueryPoint) -> Result<(), LiquidDensityError> {
        let axes = self.surface.axes();
        if point.axis_quantities() != self.context.axis_quantities()
            || point.axes().keys().ne(self.context.axes().keys())
            || point.axes().iter().any(|(axis, value)| {
                !axes.contains(axis) && value.to_bits() != self.context.axes()[axis].to_bits()
            })
        {
            return Err(invalid(
                "liquid density response cannot change fixed composition or other source context",
            ));
        }
        validate_coordinates(
            point,
            &axes[self.temperature_index],
            &axes[1 - self.temperature_index],
        )
    }
}

fn validate_coordinates(
    point: &QueryPoint,
    temperature: &str,
    pressure: &str,
) -> Result<(), LiquidDensityError> {
    let matches_axis = |axis: &str, dims: Dims, kind: QuantityKind| {
        point
            .axes()
            .get(axis)
            .is_some_and(|value| value.is_finite())
            && point.axis_quantities().get(axis).is_some_and(|quantity| {
                *quantity == QuantitySpec::dimensional(dims)
                    || *quantity
                        == QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
            })
    };
    if temperature == pressure
        || !matches_axis(
            temperature,
            Temperature::DIMS,
            QuantityKind::AbsoluteTemperature,
        )
        || !matches_axis(pressure, Pressure::DIMS, QuantityKind::Pressure)
        || point.axes()[temperature] <= 0.0
    {
        return Err(invalid(
            "liquid density response needs distinct typed absolute-temperature and pressure coordinates",
        ));
    }
    Ok(())
}

fn invalid(quantity: &'static str) -> LiquidDensityError {
    MaterialStatePointError::InvalidDerived { quantity }.into()
}

impl FluidAcoustics for ResolvedLiquidDensityState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.thermoviscous.acoustic())
    }
}
impl FluidTransport for ResolvedLiquidDensityState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.thermoviscous.transport())
    }
}
impl ThermalFluidTransport for ResolvedLiquidDensityState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        Ok(self.thermal)
    }
}
impl ThermoviscousFluid for ResolvedLiquidDensityState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        Ok(self.thermoviscous)
    }
}

#[cfg(test)]
mod tests;
