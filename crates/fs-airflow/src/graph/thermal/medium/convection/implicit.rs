//! Original liquid and gas transport on the native nonlinear exterior surface.
//!
//! The existing correlation owner supplies Nu and its analytic input partials.
//! Prepared source coefficients supply the temperature chain rule. The native
//! integrated boundary retains material/coating resistance, both bath powers,
//! bulk storage, surface elimination and time integration.

mod evaluator;
mod source;

/// Original fixed-composition gas properties on the native conducting surface.
pub mod gas;
pub use gas::ImplicitGasConvection;

/// Prescribed mass flow resolved with the original actual bulk-fluid density.
pub mod mass_flow;
pub use mass_flow::{MassFlowConvectionResponse, PrescribedMassFlowConvection};

use super::deposition::{
    ConvectionSourceError,
    source::{ConvectionPropertyTemperature, SourceCellConvectionLaw},
};
use super::*;
use fs_convection::developing::{
    FlatPlateSectionError, TurbulentTubePropertyCorrection, evaluate_circular_duct_section,
    evaluate_flat_plate_section, evaluate_turbulent_circular_duct_section,
};
use fs_matdb::PropertyLineSlice;
use fs_material::liquid::VOLUMETRIC_EXPANSION_PROPERTY;
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::MaterialStatePointError;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitGrayExteriorFluxState, ImplicitGrayExteriorLaw, ImplicitGrayExteriorState,
    ImplicitGraySurfaceState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use source::{PreparedScalar, TransportSample};

/// Immutable physical geometry, property convention and bounded source support.
#[derive(Clone, Copy, Debug)]
pub struct ExteriorConvectionSettings {
    /// Supported fluid-property temperatures, not bulk-solid temperatures.
    pub fluid_temperature_bounds_k: [f64; 2],
    /// Supported actual exterior temperatures; the boundary intersects these
    /// with its physical surface bracket. The two intervals may differ.
    pub surface_temperature_bounds_k: [f64; 2],
    pub characteristic_length: Length,
    pub property_temperature: ConvectionPropertyTemperature,
    /// Forced, finite plate/tube section, forced-sphere, natural-geometry or
    /// laminar aiding mixed law.
    pub law: SourceCellConvectionLaw,
    /// Constant total hemispherical gray emissivity on the same exterior.
    pub emissivity: f64,
}

/// Original source/domain failures remain observable before native lowering.
#[derive(Debug)]
pub enum ExteriorConvectionError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Material(MaterialError),
    Convection(ConvectionSourceError),
    Correlation(CorrelationError),
    Thermal(DepositionError),
    /// An input varies with Ts but the original card supplies no local partial.
    MissingDerivative {
        axis: &'static str,
    },
}
impl fmt::Display for ExteriorConvectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "exterior convection: {reason}"),
            Self::Source(error) => error.fmt(f),
            Self::Material(error) => error.fmt(f),
            Self::Convection(error) => error.fmt(f),
            Self::Correlation(error) => error.fmt(f),
            Self::Thermal(error) => error.fmt(f),
            Self::MissingDerivative { axis } => {
                write!(f, "exterior convection needs the original {axis} partial")
            }
        }
    }
}
impl std::error::Error for ExteriorConvectionError {}
impl From<MaterialStatePointError> for ExteriorConvectionError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<MaterialError> for ExteriorConvectionError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CorrelationError> for ExteriorConvectionError {
    fn from(error: CorrelationError) -> Self {
        Self::Correlation(error)
    }
}
impl From<FlatPlateSectionError> for ExteriorConvectionError {
    fn from(error: FlatPlateSectionError) -> Self {
        match error {
            FlatPlateSectionError::Correlation(error) => Self::Correlation(error),
            FlatPlateSectionError::UnsupportedCorrelation(_) => {
                Self::Invalid("finite plate strips require an original flat-plate mean correlation")
            }
        }
    }
}
impl From<ConvectionSourceError> for ExteriorConvectionError {
    fn from(error: ConvectionSourceError) -> Self {
        Self::Convection(error)
    }
}

enum Expansion {
    Absent,
    Density,
    Source(PreparedScalar),
}

/// One resolved coefficient and its exact local surface-temperature response.
/// Correlation evidence is retained; the tangent is not a material uncertainty
/// bound or a certification of the empirical correlation's physical accuracy.
#[derive(Debug)]
pub struct ExteriorConvectionResponse {
    pub property_temperature_k: f64,
    pub transport: FluidTransportState,
    /// The original liquid at the actual exterior temperature when the
    /// Gnielinski wall-Prandtl correction is explicitly enabled. This is a
    /// property query, never a replacement for the bulk fluid or solid storage.
    pub surface_transport: Option<FluidTransportState>,
    /// Explicit correction on the actual turbulent section, including its
    /// held bulk/wall inputs. Original cumulative mean evidence remains K=1.
    pub turbulent_property_correction: Option<TurbulentTubePropertyCorrection>,
    pub surface_dynamic_viscosity_pa_s: Option<f64>,
    /// Original correlation evaluation; for a finite section this is the
    /// downstream cumulative mean before an optional wall-property correction,
    /// not the section's actual Nusselt number.
    pub nusselt: NusseltEvaluation,
    /// Actual plate-strip geometry, upstream mean and strip-average Nu.
    pub plate_section: Option<PlateSectionDetails>,
    /// Actual circular-tube geometry, upstream mean and section-average Nu.
    pub tube_section: Option<TubeSectionDetails>,
    pub coefficient: Evidence<HeatTransferCoefficient>,
    pub coefficient_surface_derivative_w_m2_k2: f64,
}

/// Compatible name for the original liquid source/correlation error.
pub use ExteriorConvectionError as LiquidExteriorConvectionError;
/// Compatible name for the original liquid coefficient and tangent response.
pub use ExteriorConvectionResponse as LiquidExteriorConvectionResponse;
/// Compatible name for the original liquid settings. Gas and liquid sources
/// share the same physical geometry, property convention and native boundary.
pub use ExteriorConvectionSettings as LiquidExteriorConvectionSettings;

/// Source-backed liquid convection for `IntegratedCaloricGrayBoundary`.
///
/// The original liquid card, selected claims, pressure, composition and all
/// other source coordinates remain fixed. Preparation accepts constant scalars
/// with explicit thermal support and original linear-in-temperature curves.
/// It retains source receipts and native segment coefficients once; surface
/// trials perform binary searches and reuse the existing correlation evaluator,
/// without querying/cloning the material database or fitting a replacement law.
///
/// Dynamic and kinematic viscosity are both supported. Mean-film sampling
/// includes dTfilm/dTs = 1/2 in EVERY transport term. Sieder-Tate and Whitaker
/// instead retain bulk Re/Pr/k and differentiate only the same fluid's wall
/// viscosity. Natural convection uses signed sourced beta, the full Ts-Tbath
/// contrast and the existing geometry selector.
/// Aiding mixed convection retains the same source-owned expansion and the
/// caller's signed flow and density-departure bound. Its forced/natural
/// components and analytic partials remain with the original correlation.
/// Finite plate strips retain the physical leading edge, both original means
/// and the section coefficient's own Re/Pr partials. Film derivatives act on
/// that coefficient, including its conductivity. Circular sections retain bulk
/// properties (K=1 for turbulent flow by default) and zero surface derivative.
/// [`Self::with_gnielinski_wall_correction`] opts turbulent tube sections into
/// the original liquid wall-Prandtl correction and its complete source tangent.
/// Rectangular four-wall sections require a finite area owner and are not
/// pointwise laws.
///
/// A genuine source kink or unavailable active correlation partial refuses.
/// Source interval endpoints use their inward derivative. At equal sphere/bath
/// temperatures the native flux interface retains the original Nu=2 limit and
/// finite heat tangent, while the coefficient-only interface refuses singular
/// dh/dTs. Mixed-flow equilibrium similarly retains the pure-forced heat
/// tangent; at rest the heat and its tangent vanish without a correlation
/// evaluation or a conduction floor. Other natural cards retain their positive
/// minimum Rayleigh numbers.
/// No extrapolated equilibrium floor, regime blend, finite-difference Jacobian,
/// surface storage or second heat solve is introduced.
pub struct ImplicitLiquidConvection<'a, 'm> {
    caloric: &'a dyn ThermalLaw,
    energy_bounds_j_kg: [f64; 2],
    liquid: &'a LiquidTemperatureLaw<'m>,
    settings: ExteriorConvectionSettings,
    transport: [PreparedScalar; 4],
    expansion: Expansion,
    gnielinski_wall_correction: bool,
}

impl<'a, 'm> ImplicitLiquidConvection<'a, 'm> {
    /// Bind the SAME caloric object later supplied to the native boundary.
    /// Every transport source must support the declared full fluid interval;
    /// actual film/wall queries and correlation regimes are checked per trial.
    pub fn try_new(
        caloric: &'a dyn ThermalLaw,
        energy_bounds_j_kg: [f64; 2],
        liquid: &'a LiquidTemperatureLaw<'m>,
        settings: ExteriorConvectionSettings,
    ) -> Result<Self, ExteriorConvectionError> {
        use ExteriorConvectionError::Invalid;
        evaluator::admit(caloric, energy_bounds_j_kg, settings)?;
        let requirements = liquid.reference_state().requirements().properties();
        let make = |index: usize| {
            PreparedScalar::new(
                liquid,
                requirements[index].name(),
                settings.fluid_temperature_bounds_k,
                true,
            )
        };
        let transport = [make(0)?, make(1)?, make(2)?, make(3)?];
        let expansion = if let Some(original) = liquid.reference_state().density_expansion() {
            let knots = &original.source().knots;
            if settings.fluid_temperature_bounds_k[0] < knots[0].0
                || settings.fluid_temperature_bounds_k[1] > knots[knots.len() - 1].0
            {
                return Err(Invalid(
                    "fluid interval leaves the original density-expansion path",
                ));
            }
            Expansion::Density
        } else if liquid
            .reference_state()
            .material()
            .property(VOLUMETRIC_EXPANSION_PROPERTY)
            .is_some()
        {
            Expansion::Source(PreparedScalar::new(
                liquid,
                VOLUMETRIC_EXPANSION_PROPERTY,
                settings.fluid_temperature_bounds_k,
                false,
            )?)
        } else {
            Expansion::Absent
        };
        if matches!(
            settings.law,
            SourceCellConvectionLaw::VerticalPlate { .. }
                | SourceCellConvectionLaw::Natural { .. }
                | SourceCellConvectionLaw::AidingVerticalPlate { .. }
        ) && matches!(expansion, Expansion::Absent)
        {
            return Err(Invalid(
                "natural convection requires original sourced expansion",
            ));
        }
        // An augmented material bundle must not silently lose an additional
        // property's domain when replaced by a prepared transport path.
        let expected = 4 + usize::from(matches!(expansion, Expansion::Source(_)));
        if liquid.reference_state().material().properties().len() != expected {
            return Err(Invalid(
                "implicit liquid transport requires its original four-property bundle and optional expansion source",
            ));
        }
        Ok(Self {
            caloric,
            energy_bounds_j_kg,
            liquid,
            settings,
            transport,
            expansion,
            gnielinski_wall_correction: false,
        })
    }

    /// Apply Gnielinski's liquid property correction `(Pr_bulk/Pr_wall)^0.11`
    /// to the existing developing turbulent circular-tube section.
    ///
    /// The original source supplies wall Cp, conductivity and dynamic viscosity
    /// at the ACTUAL eliminated exterior temperature. Kinematic viscosity
    /// retains its original density factor. All three source derivatives enter
    /// the native surface tangent, while bulk Re, Pr, k and speed stay fixed.
    /// The same correction survives a [`PrescribedMassFlowConvection`] wrapper.
    ///
    /// This is explicit model selection, not an automatic correction for every
    /// material. The source correlation's Newtonian single-phase assumptions,
    /// property-ratio admission and section-model limitations still apply.
    /// Bulk and wall property support must both be covered by the already
    /// prepared fluid interval. Unsupported laws and source kinks refuse;
    /// no film substitution, extra heat credit or boundary solve is introduced.
    pub fn with_gnielinski_wall_correction(mut self) -> Result<Self, ExteriorConvectionError> {
        if !matches!(
            self.settings.law,
            SourceCellConvectionLaw::ForcedTubeSection { .. }
        ) || self.settings.property_temperature != ConvectionPropertyTemperature::Reservoir
        {
            return Err(ExteriorConvectionError::Invalid(
                "liquid Gnielinski wall correction requires a bulk-property turbulent circular-tube section",
            ));
        }
        let [lo, hi] = self.settings.fluid_temperature_bounds_k;
        let [wall_lo, wall_hi] = self.settings.surface_temperature_bounds_k;
        if wall_lo < lo || wall_hi > hi {
            return Err(ExteriorConvectionError::Invalid(
                "wall-Prandtl correction requires original liquid support over the whole surface interval",
            ));
        }
        self.gnielinski_wall_correction = true;
        Ok(self)
    }

    /// Immutable original liquid owner, including its card and complete context.
    pub const fn liquid(&self) -> &LiquidTemperatureLaw<'m> {
        self.liquid
    }

    /// Original density/Cp/k/viscosity source lines, in that physical order.
    pub fn transport_sources(&self) -> impl Iterator<Item = &PropertyLineSlice> {
        self.transport.iter().map(|source| &source.source)
    }

    /// Independently sourced beta, when selected instead of the density tangent.
    pub fn expansion_source(&self) -> Option<&PropertyLineSlice> {
        match &self.expansion {
            Expansion::Source(source) => Some(&source.source),
            _ => None,
        }
    }

    /// Evaluate physical Ts and bath T with typed source/correlation failures.
    /// The historical native `gas_temperature_k` label also denotes a liquid bath.
    pub fn evaluate(
        &self,
        surface_temperature_k: f64,
        bath_temperature_k: f64,
    ) -> Result<ExteriorConvectionResponse, ExteriorConvectionError> {
        evaluator::evaluate(self, surface_temperature_k, bath_temperature_k)
    }
}

impl evaluator::ExteriorFluid for ImplicitLiquidConvection<'_, '_> {
    fn caloric(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn energy_bounds(&self) -> [f64; 2] {
        self.energy_bounds_j_kg
    }
    fn settings(&self) -> ExteriorConvectionSettings {
        self.settings
    }

    fn turbulent_wall_correction(
        &self,
        temperature: f64,
    ) -> Result<Option<evaluator::TurbulentWallCorrection>, ExteriorConvectionError> {
        if self.gnielinski_wall_correction {
            self.sample_transport(temperature, true)
                .map(|sample| Some(evaluator::TurbulentWallCorrection::Liquid(sample)))
        } else {
            Ok(None)
        }
    }

    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, ExteriorConvectionError> {
        let density_derivative = derivatives || matches!(self.expansion, Expansion::Density);
        let (rho, drho) = self.transport[0].sample(temperature, density_derivative)?;
        let (cp, dcp) = self.transport[1].sample(temperature, derivatives)?;
        let (k, dk) = self.transport[2].sample(temperature, derivatives)?;
        let (viscosity, slope) = self.transport[3].sample(temperature, derivatives)?;
        let (mu, dmu) = match self.liquid.viscosity_basis() {
            LiquidViscosityBasis::Dynamic => (viscosity, slope),
            LiquidViscosityBasis::Kinematic => (
                finite(rho * viscosity)?,
                finite(drho * viscosity + rho * slope)?,
            ),
        };
        Ok(TransportSample {
            state: FluidTransportState::try_new(rho, cp, k, mu)?,
            derivatives: [drho, dcp, dk, dmu],
        })
    }

    fn wall_viscosity(
        &self,
        temperature: f64,
    ) -> Result<(f64, f64, bool), ExteriorConvectionError> {
        let (viscosity, derivative) = self.transport[3].sample(temperature, true)?;
        match self.liquid.viscosity_basis() {
            LiquidViscosityBasis::Dynamic => Ok((viscosity, derivative, derivative != 0.0)),
            LiquidViscosityBasis::Kinematic => {
                let (rho, drho) = self.transport[0].sample(temperature, true)?;
                Ok((
                    finite(rho * viscosity)?,
                    finite(drho * viscosity + rho * derivative)?,
                    drho != 0.0 || derivative != 0.0,
                ))
            }
        }
    }

    fn expansion(
        &self,
        temperature: f64,
        sample: &TransportSample,
        derivative: bool,
    ) -> Result<(f64, f64), ExteriorConvectionError> {
        match &self.expansion {
            Expansion::Absent => Err(ExteriorConvectionError::Invalid(
                "missing sourced expansion",
            )),
            Expansion::Source(source) => source.sample(temperature, derivative),
            Expansion::Density => {
                let beta = finite(-sample.derivatives[0] / sample.state.density_kg_m3())?;
                // rho is affine on the original segment, so rho''=0.
                Ok((beta, finite(beta * beta)?))
            }
        }
    }
}

impl ImplicitGrayExteriorLaw for ImplicitLiquidConvection<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        self.energy_bounds_j_kg
    }
    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        self.settings.surface_temperature_bounds_k
    }
    fn flux_state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorFluxState, DepositionError> {
        evaluator::flux_state_at_surface(self, energy, surface, bath)
    }
    fn state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        evaluator::state_at_surface(self, energy, surface, bath)
    }
}

fn finite(value: f64) -> Result<f64, ExteriorConvectionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ExteriorConvectionError::Invalid(
            "nonfinite source or tangent arithmetic",
        ))
    }
}
fn positive_interval([lo, hi]: [f64; 2]) -> bool {
    lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo
}
