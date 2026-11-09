//! An isothermal bare bank exchanging heat with a finite liquid stream.
//!
//! The original Cp source owns the entire heat integral and its inverse. A
//! bounded mean-temperature solve refreshes the SAME liquid's rho, Cp, k and
//! viscosity at (T_in + T_out)/2; its surface Pr belongs at T_wall. The
//! representative approach speed is mdot/(rho_mean A_frontal), conserving the
//! same mass flux as the inlet. Neither velocity uses heated surface area.
//! The effective Robin reference comes from the native caloric exchange mean,
//! which generally differs from the arithmetic property temperature.
//!
//! The bank is steady, uniformly isothermal and free of bypass. Its long bare
//! tubes expose lateral area only; tube ends, fins, wall/fluid storage, pressure
//! losses, natural convection and a second finite stream are not modelled.

mod mean;

use super::{MaterialTubeBankError, evaluate_tube_bank_from_material};
use fs_convection::tube_bank::{TubeBankError, TubeBankEvaluation, TubeBankGeometry};
use fs_convection::{CorrelationRobinBoundary, ThermalConductivity};
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, InterpolationPolicy, PropertyValue};
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::state_point::MaterialStatePointError;
use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::isobaric::{IsobaricHeatCapacityLaw, LinearHeatExchange};
use fs_qty::{Length, Temperature};

/// Lateral tube area and frontal area from one original bank geometry.
///
/// Each row has `tubes_per_row` full pitch cells across its width. The frontal
/// area is N_transverse * S_T * span, with half-pitch lateral margins or a
/// transverse repeating cell. Clearance/bypass and exposed tube ends cannot be
/// added to this area. An axial span is not a fin height or a short-pin model.
/// For staggered rows this is a transverse repeating strip; shifted edge tubes
/// are not a separately resolved finite-bank envelope or side-bypass path.
/// The caller must establish negligible tube-end influence; this area identity
/// does not establish a quantitative aspect-ratio validity range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeBankSurface {
    bank: TubeBankGeometry,
    tubes_per_row: u32,
    span: Length,
    frontal_area_m2: f64,
    surface_area_m2: f64,
}

impl TubeBankSurface {
    pub fn new(
        bank: TubeBankGeometry,
        tubes_per_row: u32,
        span: Length,
    ) -> Result<Self, TubeBankStreamError> {
        positive(span.value(), "tube span must be positive")?;
        if tubes_per_row == 0 {
            return Err(TubeBankStreamError::Invalid(
                "a bank needs tubes in every row",
            ));
        }
        let length = positive(
            f64::from(tubes_per_row) * span.value(),
            "tube length overflow",
        )?;
        let frontal_area_m2 = positive(length * bank.transverse_pitch().value(), "frontal area")?;
        let surface_area_m2 = positive(
            length * (f64::from(bank.rows()) * core::f64::consts::PI * bank.diameter().value()),
            "lateral tube area",
        )?;
        Ok(Self {
            bank,
            tubes_per_row,
            span,
            frontal_area_m2,
            surface_area_m2,
        })
    }
    pub const fn bank(self) -> TubeBankGeometry {
        self.bank
    }
    pub const fn tubes_per_row(self) -> u32 {
        self.tubes_per_row
    }
    pub const fn span(self) -> Length {
        self.span
    }
    pub const fn frontal_area_m2(self) -> f64 {
        self.frontal_area_m2
    }
    pub const fn surface_area_m2(self) -> f64 {
        self.surface_area_m2
    }
}

/// Bounds property iteration, source traversal and the independent watt audits.
#[derive(Clone, Copy, Debug)]
pub struct TubeBankStreamConfig {
    pub temperature_tolerance_k: f64,
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
    pub maximum_iterations: usize,
    pub maximum_source_knots: usize,
}
impl Default for TubeBankStreamConfig {
    fn default() -> Self {
        Self {
            temperature_tolerance_k: 1e-9,
            absolute_heat_tolerance_w: 1e-7,
            relative_heat_tolerance: 1e-10,
            maximum_iterations: 64,
            maximum_source_knots: 16_384,
        }
    }
}

#[derive(Debug)]
pub enum TubeBankStreamError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Caloric(ThermalTransportError),
    Convection(MaterialTubeBankError),
    HeatBalance { residual_w: f64, tolerance_w: f64 },
    NonConvergence { iterations: usize, residual_k: f64 },
    Cancelled,
}
impl core::fmt::Display for TubeBankStreamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite tube-bank stream: {self:?}")
    }
}
impl std::error::Error for TubeBankStreamError {}
impl From<MaterialStatePointError> for TubeBankStreamError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Source(value)
    }
}
impl From<ThermalTransportError> for TubeBankStreamError {
    fn from(value: ThermalTransportError) -> Self {
        Self::Caloric(value)
    }
}
impl From<MaterialTubeBankError> for TubeBankStreamError {
    fn from(value: MaterialTubeBankError) -> Self {
        Self::Convection(value)
    }
}

/// One immutable source/geometry/flow binding. Solves publish only complete states.
pub struct LiquidTubeBankStream<'a, 'card> {
    liquid: &'a LiquidTemperatureLaw<'card>,
    caloric: &'a IsobaricHeatCapacityLaw,
    surface: TubeBankSurface,
    inlet_temperature_k: f64,
    mass_flow_kg_s: f64,
}

/// Original material states and one reciprocal heat rate at the accepted endpoint.
#[derive(Debug)]
pub struct TubeBankStreamState<'a> {
    caloric: &'a IsobaricHeatCapacityLaw,
    inlet_source: ResolvedLiquidState,
    bulk_source: ResolvedLiquidState,
    surface_source: ResolvedLiquidState,
    convection: TubeBankEvaluation,
    surface: TubeBankSurface,
    inlet_temperature_k: f64,
    wall_temperature_k: f64,
    property_temperature_k: f64,
    mass_flow_kg_s: f64,
    /// Actual inlet mass-flow velocity before the bank's blockage correction.
    inlet_speed_m_s: f64,
    /// Representative unblocked velocity at the mean bulk source density.
    approach_speed_m_s: f64,
    conductance_w_k: f64,
    /// Original caloric inverse, relative enthalpy and effective Robin mean.
    exchange: LinearHeatExchange,
    /// Positive from the wall into the stream. The wall loses this SAME heat.
    heat_into_fluid_w: f64,
    /// Original Cp integral at represented T_out minus retained relative heat.
    endpoint_energy_residual_w: f64,
    robin_energy_residual_w: f64,
    property_residual_k: f64,
    iterations: usize,
}
impl<'a> TubeBankStreamState<'a> {
    /// Original integrated Cp law, including its source span and receipts.
    pub const fn caloric(&self) -> &'a IsobaricHeatCapacityLaw {
        self.caloric
    }
    pub const fn inlet_source(&self) -> &ResolvedLiquidState {
        &self.inlet_source
    }
    pub const fn bulk_source(&self) -> &ResolvedLiquidState {
        &self.bulk_source
    }
    pub const fn surface_source(&self) -> &ResolvedLiquidState {
        &self.surface_source
    }
    pub const fn convection(&self) -> &TubeBankEvaluation {
        &self.convection
    }
    pub const fn surface(&self) -> TubeBankSurface {
        self.surface
    }
    pub const fn inlet_temperature_k(&self) -> f64 {
        self.inlet_temperature_k
    }
    pub const fn wall_temperature_k(&self) -> f64 {
        self.wall_temperature_k
    }
    /// Arithmetic inlet/outlet mean used to evaluate the bank's material properties.
    pub const fn property_temperature_k(&self) -> f64 {
        self.property_temperature_k
    }
    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.mass_flow_kg_s
    }
    pub const fn inlet_speed_m_s(&self) -> f64 {
        self.inlet_speed_m_s
    }
    pub const fn approach_speed_m_s(&self) -> f64 {
        self.approach_speed_m_s
    }
    pub const fn conductance_w_k(&self) -> f64 {
        self.conductance_w_k
    }
    pub const fn exchange(&self) -> &LinearHeatExchange {
        &self.exchange
    }
    pub const fn heat_into_fluid_w(&self) -> f64 {
        self.heat_into_fluid_w
    }
    pub const fn endpoint_energy_residual_w(&self) -> f64 {
        self.endpoint_energy_residual_w
    }
    pub const fn robin_energy_residual_w(&self) -> f64 {
        self.robin_energy_residual_w
    }
    pub const fn property_residual_k(&self) -> f64 {
        self.property_residual_k
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub fn heat_into_wall_w(&self) -> f64 {
        -self.heat_into_fluid_w
    }

    /// Actual coefficient and caloric mean; no caller-substituted bulk reference.
    pub fn robin_boundary(&self) -> Result<CorrelationRobinBoundary, TubeBankError> {
        self.convection.robin_boundary(
            ThermalConductivity::new(self.bulk_source.thermal_conductivity_w_m_k()),
            Temperature::new(self.exchange.mean_temperature_k),
        )
    }
}

impl<'a, 'card> LiquidTubeBankStream<'a, 'card> {
    pub fn new(
        liquid: &'a LiquidTemperatureLaw<'card>,
        caloric: &'a IsobaricHeatCapacityLaw,
        surface: TubeBankSurface,
        inlet_temperature_k: f64,
        mass_flow_kg_s: f64,
    ) -> Result<Self, TubeBankStreamError> {
        positive(inlet_temperature_k, "inlet temperature")?;
        // A stagnant bank needs a separate natural/conductive heat owner; its
        // forced-throughflow correlation cannot supply a zero-speed limit.
        positive(mass_flow_kg_s, "forced bank requires positive throughflow")?;
        let reference = liquid.reference_point();
        let endpoint = &caloric.query_endpoints()[0];
        if caloric.material_card_identity() != liquid.card().content_hash()
            || caloric.temperature_axis() != liquid.temperature_axis()
            || endpoint.axis_quantities() != reference.axis_quantities()
            || !endpoint
                .axes()
                .iter()
                .filter(|(k, _)| k.as_str() != liquid.temperature_axis())
                .eq(reference
                    .axes()
                    .iter()
                    .filter(|(k, _)| k.as_str() != liquid.temperature_axis()))
        {
            return Err(TubeBankStreamError::Invalid(
                "calorics and transport need the same original card and fixed source context",
            ));
        }
        caloric.specific_heat_cp_j_kg_k(inlet_temperature_k)?;
        Ok(Self {
            liquid,
            caloric,
            surface,
            inlet_temperature_k,
            mass_flow_kg_s,
        })
    }

    pub const fn liquid(&self) -> &'a LiquidTemperatureLaw<'card> {
        self.liquid
    }
    pub const fn caloric(&self) -> &'a IsobaricHeatCapacityLaw {
        self.caloric
    }
    pub const fn surface(&self) -> TubeBankSurface {
        self.surface
    }

    /// Solve one uniformly isothermal bank, refreshing original liquid sources.
    ///
    /// Every source must be continuous over the inlet-to-wall interval. The
    /// original Nu domain must hold at every attempted mean-property point;
    /// refusal never clips Re/Pr or substitutes a previous coefficient. The
    /// reported temperature residual is a numerical closure tolerance, not an
    /// empirical accuracy or uniqueness claim for the mean-property model.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        wall_temperature_k: f64,
        config: TubeBankStreamConfig,
    ) -> Result<TubeBankStreamState<'a>, TubeBankStreamError> {
        poll(cx)?;
        positive(wall_temperature_k, "wall temperature")?;
        positive(config.temperature_tolerance_k, "temperature tolerance")?;
        if !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
            || config.maximum_iterations == 0
            || config.maximum_source_knots < 2
        {
            return Err(TubeBankStreamError::Invalid(
                "invalid stream work or heat-balance bounds",
            ));
        }
        self.admit_sources(cx, wall_temperature_k, config.maximum_source_knots)?;
        let inlet = self
            .liquid
            .at_temperature(Temperature::new(self.inlet_temperature_k))?;
        let wall = self
            .liquid
            .at_temperature(Temperature::new(wall_temperature_k))?;
        let mass_flux = positive(
            self.mass_flow_kg_s / self.surface.frontal_area_m2,
            "frontal mass flux",
        )?;
        let inlet_speed = positive(
            mass_flux / inlet.density_kg_m3(),
            "inlet mass-flow velocity",
        )?;
        let solved = mean::solve(
            self.inlet_temperature_k,
            wall_temperature_k,
            config.temperature_tolerance_k,
            config.maximum_iterations,
            || poll(cx),
            |mean| {
                let bulk = self.liquid.at_temperature(Temperature::new(mean))?;
                let speed = positive(mass_flux / bulk.density_kg_m3(), "mean mass-flow velocity")?;
                let convection =
                    evaluate_tube_bank_from_material(&bulk, &wall, speed, self.surface.bank)?;
                let conductance = positive(
                    convection.coefficient().value.value() * self.surface.surface_area_m2,
                    "bank conductance",
                )?;
                let exposure = positive(conductance / self.mass_flow_kg_s, "bank exposure")?;
                let exchange = self.caloric.linear_heat_exchange(
                    self.inlet_temperature_k,
                    wall_temperature_k,
                    exposure,
                )?;
                let outlet = exchange.temperature_k;
                let evaluation = convection.evaluation;
                Ok(((bulk, evaluation, speed, conductance, exchange), outlet))
            },
        )
        .map_err(|error| match error {
            mean::MeanError::Callback(error) => error,
            mean::MeanError::InvalidOutlet => TubeBankStreamError::Invalid(
                "source exchange leaves the passive temperature bracket",
            ),
            mean::MeanError::Unresolved {
                iterations,
                residual,
            } => TubeBankStreamError::NonConvergence {
                iterations,
                residual_k: residual,
            },
        })?;
        let (bulk, evaluation, speed, conductance, exchange) = solved.value;
        let heat = self.mass_flow_kg_s * exchange.specific_enthalpy_change_j_kg;
        if exchange.specific_enthalpy_change_j_kg != 0.0 && heat == 0.0 {
            return Err(TubeBankStreamError::Invalid(
                "nonzero bank heat is unrepresentable",
            ));
        }
        let endpoint_heat = self.mass_flow_kg_s
            * self
                .caloric
                .specific_enthalpy_change_j_kg(self.inlet_temperature_k, exchange.temperature_k)?;
        let robin_heat = conductance * (wall_temperature_k - exchange.mean_temperature_k);
        let tolerance =
            config.absolute_heat_tolerance_w + config.relative_heat_tolerance * heat.abs();
        let endpoint_residual = endpoint_heat - heat;
        let robin_residual = robin_heat - heat;
        for error in [endpoint_residual, robin_residual] {
            if !heat.is_finite()
                || !tolerance.is_finite()
                || !error.is_finite()
                || error.abs() > tolerance
            {
                return Err(TubeBankStreamError::HeatBalance {
                    residual_w: error,
                    tolerance_w: tolerance,
                });
            }
        }
        poll(cx)?;
        Ok(TubeBankStreamState {
            caloric: self.caloric,
            inlet_source: inlet,
            bulk_source: bulk,
            surface_source: wall,
            convection: evaluation,
            surface: self.surface,
            inlet_temperature_k: self.inlet_temperature_k,
            wall_temperature_k,
            property_temperature_k: solved.mean,
            mass_flow_kg_s: self.mass_flow_kg_s,
            inlet_speed_m_s: inlet_speed,
            approach_speed_m_s: speed,
            conductance_w_k: conductance,
            exchange,
            heat_into_fluid_w: heat,
            endpoint_energy_residual_w: endpoint_residual,
            robin_energy_residual_w: robin_residual,
            property_residual_k: solved.residual,
            iterations: solved.iterations,
        })
    }

    fn admit_sources(
        &self,
        cx: &Cx<'_>,
        wall: f64,
        limit: usize,
    ) -> Result<(), TubeBankStreamError> {
        // Wall-fluid Pr needs the same Cp claim at the surface. Requiring the
        // complete caloric bracket also prevents trials from escaping a shorter
        // prepared heat-capacity path even when the source card is wider.
        self.caloric.specific_heat_cp_j_kg_k(wall)?;
        let original = self.liquid.reference_state();
        let properties = original.requirements().properties();
        let cp = original.material().property(properties[1].name()).ok_or(
            TubeBankStreamError::Invalid("original liquid Cp claim is absent"),
        )?;
        if self.caloric.source_line().knots.len() > limit {
            return Err(TubeBankStreamError::Invalid(
                "caloric source knot budget exceeded",
            ));
        }
        for (i, (_, answer)) in self.caloric.source_line().knots.iter().enumerate() {
            if i % 64 == 0 {
                poll(cx)?;
            }
            if answer.receipt.selected != cp.answer().receipt.selected {
                return Err(TubeBankStreamError::Invalid(
                    "transport and calorics must select the same original Cp claim",
                ));
            }
        }
        for property in properties {
            poll(cx)?;
            let selected = original
                .material()
                .property(property.name())
                .ok_or(TubeBankStreamError::Invalid(
                    "original liquid transport claim is absent",
                ))?
                .answer()
                .receipt
                .selected;
            let claim =
                self.liquid
                    .card()
                    .claims()
                    .claim(selected)
                    .ok_or(TubeBankStreamError::Invalid(
                        "original liquid claim is absent",
                    ))?;
            match &claim.value {
                PropertyValue::Scalar { .. }
                    if claim.interpolation == InterpolationPolicy::ConstantWithinValidity => {}
                PropertyValue::Curve {
                    abscissa, knots, ..
                } if abscissa == self.liquid.temperature_axis()
                    && claim.interpolation == InterpolationPolicy::LinearInside
                    && knots.len() <= limit => {}
                _ => {
                    return Err(TubeBankStreamError::Invalid(
                        "bank iteration needs bounded continuous original scalar/temperature-curve transport",
                    ));
                }
            }
            if wall != self.inlet_temperature_k {
                let point = |temperature| {
                    self.liquid
                        .reference_point()
                        .clone()
                        .with(self.liquid.temperature_axis(), temperature)
                        .map_err(|source| MaterialStatePointError::Query {
                            property: property.name().to_owned(),
                            source,
                        })
                };
                let line = self
                    .liquid
                    .card()
                    .claims()
                    .query_line_typed(
                        &claim.key,
                        self.liquid.temperature_axis(),
                        &point(wall.min(self.inlet_temperature_k))?,
                        &point(wall.max(self.inlet_temperature_k))?,
                        ClaimSelection::Pinned(selected),
                    )
                    .map_err(|source| MaterialStatePointError::Support {
                        property: property.name().to_owned(),
                        source: Box::new(source),
                    })?;
                if line.knots.len() > limit {
                    return Err(TubeBankStreamError::Invalid(
                        "transport source knot budget exceeded",
                    ));
                }
                for (index, (_, answer)) in line.knots.iter().enumerate() {
                    if index % 64 == 0 {
                        poll(cx)?;
                    }
                    positive(
                        answer.evidence.value.value,
                        "transport source must stay positive through the stream",
                    )?;
                }
            }
        }
        Ok(())
    }
}

fn positive(value: f64, name: &'static str) -> Result<f64, TubeBankStreamError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(TubeBankStreamError::Invalid(name))
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), TubeBankStreamError> {
    cx.checkpoint().map_err(|_| TubeBankStreamError::Cancelled)
}

#[cfg(test)]
mod tests;
