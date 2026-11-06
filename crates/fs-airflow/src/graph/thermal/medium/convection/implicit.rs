//! Original liquid transport curves on the native nonlinear exterior surface.
//!
//! The existing correlation owner supplies Nu and its analytic input partials.
//! Prepared source coefficients supply the temperature chain rule. The native
//! integrated boundary retains material/coating resistance, both bath powers,
//! bulk storage, surface elimination and time integration.

mod source;

use super::deposition::{
    ConvectionSourceError,
    source::{ConvectionPropertyTemperature, SourceCellConvectionLaw},
};
use super::*;
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
pub struct LiquidExteriorConvectionSettings {
    /// Supported fluid-property temperatures, not bulk-solid temperatures.
    pub fluid_temperature_bounds_k: [f64; 2],
    /// Supported actual exterior temperatures; the boundary intersects these
    /// with its physical surface bracket. The two intervals may differ.
    pub surface_temperature_bounds_k: [f64; 2],
    pub characteristic_length: Length,
    pub property_temperature: ConvectionPropertyTemperature,
    /// Forced, forced-sphere, natural-geometry or laminar aiding mixed law.
    /// Finite plate strips are not admitted here.
    pub law: SourceCellConvectionLaw,
    /// Constant total hemispherical gray emissivity on the same exterior.
    pub emissivity: f64,
}

/// Original source/domain failures remain observable before native lowering.
#[derive(Debug)]
pub enum LiquidExteriorConvectionError {
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
impl fmt::Display for LiquidExteriorConvectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "liquid exterior convection: {reason}"),
            Self::Source(error) => error.fmt(f),
            Self::Material(error) => error.fmt(f),
            Self::Convection(error) => error.fmt(f),
            Self::Correlation(error) => error.fmt(f),
            Self::Thermal(error) => error.fmt(f),
            Self::MissingDerivative { axis } => {
                write!(
                    f,
                    "liquid exterior convection needs the original {axis} partial"
                )
            }
        }
    }
}
impl std::error::Error for LiquidExteriorConvectionError {}
impl From<MaterialStatePointError> for LiquidExteriorConvectionError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<MaterialError> for LiquidExteriorConvectionError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CorrelationError> for LiquidExteriorConvectionError {
    fn from(error: CorrelationError) -> Self {
        Self::Correlation(error)
    }
}
impl From<ConvectionSourceError> for LiquidExteriorConvectionError {
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
pub struct LiquidExteriorConvectionResponse {
    pub property_temperature_k: f64,
    pub transport: FluidTransportState,
    pub surface_dynamic_viscosity_pa_s: Option<f64>,
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
    pub coefficient_surface_derivative_w_m2_k2: f64,
}

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
    settings: LiquidExteriorConvectionSettings,
    transport: [PreparedScalar; 4],
    expansion: Expansion,
}

impl<'a, 'm> ImplicitLiquidConvection<'a, 'm> {
    /// Bind the SAME caloric object later supplied to the native boundary.
    /// Every transport source must support the declared full fluid interval;
    /// actual film/wall queries and correlation regimes are checked per trial.
    pub fn try_new(
        caloric: &'a dyn ThermalLaw,
        energy_bounds_j_kg: [f64; 2],
        liquid: &'a LiquidTemperatureLaw<'m>,
        settings: LiquidExteriorConvectionSettings,
    ) -> Result<Self, LiquidExteriorConvectionError> {
        use LiquidExteriorConvectionError::Invalid;
        let [lo, hi] = energy_bounds_j_kg;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || !positive_interval(settings.fluid_temperature_bounds_k)
            || !positive_interval(settings.surface_temperature_bounds_k)
            || !settings.characteristic_length.value().is_finite()
            || settings.characteristic_length.value() <= 0.0
            || !settings.emissivity.is_finite()
            || !(0.0..=1.0).contains(&settings.emissivity)
        {
            return Err(Invalid(
                "finite energy and positive temperature/length bounds and passive emissivity are required",
            ));
        }
        for energy in energy_bounds_j_kg {
            let temperature = caloric
                .temperature_k(energy)
                .map_err(LiquidExteriorConvectionError::Thermal)?;
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(Invalid("the original caloric endpoints must be physical"));
            }
        }
        if matches!(
            settings.law,
            SourceCellConvectionLaw::ForcedPlateSection { .. }
        ) {
            return Err(Invalid(
                "implicit liquid exterior does not support finite plate strips",
            ));
        }
        if let SourceCellConvectionLaw::AidingVerticalPlate { settings } = settings.law {
            settings.validate()?;
        }
        settings
            .law
            .admit_geometry(settings.characteristic_length)?;
        settings
            .law
            .admit_property_temperature(settings.property_temperature)?;
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
        })
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
    ) -> Result<LiquidExteriorConvectionResponse, LiquidExteriorConvectionError> {
        use LiquidExteriorConvectionError::Invalid;
        if !surface_temperature_k.is_finite()
            || !(self.settings.surface_temperature_bounds_k[0]
                ..=self.settings.surface_temperature_bounds_k[1])
                .contains(&surface_temperature_k)
            || !bath_temperature_k.is_finite()
            || bath_temperature_k <= 0.0
        {
            return Err(Invalid(
                "exterior or bath temperature leaves its physical support",
            ));
        }
        let (property_temperature_k, film_factor) = match self.settings.property_temperature {
            ConvectionPropertyTemperature::Reservoir => (bath_temperature_k, 0.0),
            ConvectionPropertyTemperature::FilmMean => {
                (0.5 * surface_temperature_k + 0.5 * bath_temperature_k, 0.5)
            }
        };
        let sample = self.sample_transport(property_temperature_k, film_factor != 0.0)?;
        let state = sample.state;
        let length = self.settings.characteristic_length;
        let rho = state.density_kg_m3();
        let cp = state.specific_heat_j_kg_k();
        let k = state.thermal_conductivity_w_m_k();
        let mu = state.dynamic_viscosity_pa_s();
        let [drho, dcp, dk, dmu] = sample.derivatives.map(|slope| slope * film_factor);
        let dpr = finite(state.prandtl() * (dmu / mu + dcp / cp - dk / k))?;
        let varying_transport =
            film_factor != 0.0 && sample.derivatives.iter().any(|slope| *slope != 0.0);
        let mut wall_viscosity = None;
        let mut varying_wall_viscosity = false;
        let (nusselt, dnu) = match self.settings.law {
            SourceCellConvectionLaw::Forced {
                correlation,
                speed_m_s,
                settings,
            } => {
                if correlation == CorrelationId::DittusBoelter
                    && ((surface_temperature_k < bath_temperature_k
                        && settings.direction == ThermalDirection::HeatingFluid)
                        || (surface_temperature_k > bath_temperature_k
                            && settings.direction == ThermalDirection::CoolingFluid))
                {
                    return Err(Invalid(
                        "Dittus-Boelter direction disagrees with actual surface/bath heat flow",
                    ));
                }
                let re = state.reynolds(speed_m_s, length.value())?;
                let dre = finite(re * (drho / rho - dmu / mu))?;
                let mut inputs = CorrelationInputs::forced(re, state.prandtl())
                    .with_direction(settings.direction);
                if let Some(ratio) = settings.length_over_hydraulic_diameter {
                    inputs = inputs.with_length_ratio(ratio);
                }
                if let Some(ratio) = settings.aspect_ratio {
                    inputs = inputs.with_aspect_ratio(ratio);
                }
                let dratio = if correlation == CorrelationId::SiederTateTurbulent {
                    let (wall, dwall, varies) = self.wall_viscosity(surface_temperature_k)?;
                    varying_wall_viscosity = varies;
                    wall_viscosity = Some(wall);
                    let ratio = finite(mu / wall)?;
                    inputs = inputs.with_viscosity_ratio(ratio);
                    finite(-ratio * dwall / wall)?
                } else {
                    0.0
                };
                let nu = inputs.evaluate_with_sensitivities(correlation)?;
                let derivative = partial(&nu, "Re", dre)?
                    + partial(&nu, "Pr", dpr)?
                    + partial(&nu, "mu_bulk_over_surface", dratio)?;
                (nu, finite(derivative)?)
            }
            SourceCellConvectionLaw::ForcedSphere { speed_m_s } => {
                let re = state.reynolds(speed_m_s, length.value())?;
                let (wall, dwall, varies) = self.wall_viscosity(surface_temperature_k)?;
                varying_wall_viscosity = varies;
                wall_viscosity = Some(wall);
                let ratio = finite(mu / wall)?;
                let nu = CorrelationInputs::forced(re, state.prandtl())
                    .with_viscosity_ratio(ratio)
                    .evaluate_with_sensitivities(CorrelationId::WhitakerSphere)?;
                let derivative =
                    partial(&nu, "mu_bulk_over_surface", finite(-ratio * dwall / wall)?)?;
                (nu, derivative)
            }
            SourceCellConvectionLaw::AidingVerticalPlate { settings } => {
                let (beta, dbeta_dt) =
                    self.expansion(property_temperature_k, &sample, film_factor != 0.0)?;
                let delta = surface_temperature_k - bath_temperature_k;
                let groups = settings.groups(
                    state,
                    beta,
                    Temperature::new(surface_temperature_k),
                    Temperature::new(bath_temperature_k),
                    length,
                )?;
                let constant_zero_drive = settings.gravity_m_s2 == 0.0
                    || (beta == 0.0 && dbeta_dt == 0.0);
                if groups.rayleigh == 0.0 && !constant_zero_drive {
                    return Err(Invalid(
                        "zero-Rayleigh mixed convection needs the native heat-tangent interface",
                    ));
                }
                let dre = finite(groups.reynolds * (drho / rho - dmu / mu))?;
                let factor = finite(
                    settings.gravity_m_s2 * length.value().powi(3)
                        / (state.kinematic_viscosity_m2_s() * state.thermal_diffusivity_m2_s()),
                )?;
                // Differentiate the signed drive before taking its magnitude;
                // neither beta nor Delta T is used as a divisor.
                let dra = finite(
                    factor
                        * (beta
                            + delta
                                * (film_factor * dbeta_dt
                                    + beta * (2.0 * drho / rho + dcp / cp - dmu / mu - dk / k)))
                        * beta.signum()
                        * delta.signum(),
                )?;
                let nu = CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl)
                    .evaluate_with_sensitivities(CorrelationId::ManningQureshiAidingVerticalPlate)?;
                require_partial(&nu, "Ra", groups.rayleigh != 0.0)?;
                let derivative = partial(&nu, "Re", dre)?
                    + partial(&nu, "Ra", dra)?
                    + partial(&nu, "Pr", dpr)?;
                (nu, finite(derivative)?)
            }
            SourceCellConvectionLaw::VerticalPlate { gravity_m_s2 }
            | SourceCellConvectionLaw::Natural { gravity_m_s2, .. } => {
                let geometry = match self.settings.law {
                    SourceCellConvectionLaw::Natural { geometry, .. } => geometry,
                    _ => NaturalConvectionGeometry::VerticalPlate,
                };
                let (beta, dbeta_dt) =
                    self.expansion(property_temperature_k, &sample, film_factor != 0.0)?;
                let delta = surface_temperature_k - bath_temperature_k;
                let ra = state.rayleigh(gravity_m_s2, beta, delta, length.value())?;
                let constant_zero_drive = geometry == NaturalConvectionGeometry::Sphere
                    && (gravity_m_s2 == 0.0
                        || (beta == 0.0 && (film_factor == 0.0 || dbeta_dt == 0.0)));
                if ra == 0.0 && !constant_zero_drive {
                    return Err(Invalid(
                        "zero-Rayleigh natural convection has no finite coefficient tangent on this interface",
                    ));
                }
                // No division by beta or Delta T: signed zero/negative source
                // expansion must not turn into a manufactured ideal-gas law.
                let factor = finite(
                    gravity_m_s2 * length.value().powi(3)
                        / (state.kinematic_viscosity_m2_s() * state.thermal_diffusivity_m2_s()),
                )?;
                let dra = finite(
                    factor
                        * (beta
                            + delta
                                * (film_factor * dbeta_dt
                                    + beta * (2.0 * drho / rho + dcp / cp - dmu / mu - dk / k))),
                )?;
                let nu = geometry.evaluate_with_sensitivities(ra, state.prandtl())?;
                // A zero first-order chain factor does not justify selecting a
                // derivative at a discontinuous natural-regime boundary.
                require_partial(&nu, "Ra", ra != 0.0)?;
                let derivative = partial(&nu, "Ra", dra * ra.signum())? + partial(&nu, "Pr", dpr)?;
                (nu, finite(derivative)?)
            }
            _ => unreachable!("constructor excludes unsupported geometry"),
        };
        // Moving source properties can have a zero aggregate group derivative
        // at a local extremum. That does not establish differentiability of a
        // missing empirical partial at a regime boundary.
        if varying_transport {
            for axis in ["Re", "Pr"] {
                // A stopped external stream stays at Re=0 regardless of the
                // fluid-property response; it needs no Reynolds derivative.
                let varies = nusselt.groups().get(axis)
                    .is_some_and(|&value| axis != "Re" || value != 0.0);
                require_partial(&nusselt, axis, varies)?;
            }
        }
        require_partial(&nusselt, "mu_bulk_over_surface", varying_wall_viscosity)?;
        let coefficient = nusselt.heat_transfer_coefficient(ThermalConductivity::new(k), length)?;
        let derivative = finite((dnu * k + nusselt.evidence().value * dk) / length.value())?;
        Ok(LiquidExteriorConvectionResponse {
            property_temperature_k,
            transport: state,
            surface_dynamic_viscosity_pa_s: wall_viscosity,
            nusselt,
            coefficient,
            coefficient_surface_derivative_w_m2_k2: derivative,
        })
    }

    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, LiquidExteriorConvectionError> {
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
    ) -> Result<(f64, f64, bool), LiquidExteriorConvectionError> {
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
    ) -> Result<(f64, f64), LiquidExteriorConvectionError> {
        match &self.expansion {
            Expansion::Absent => Err(LiquidExteriorConvectionError::Invalid(
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
        // At exact sphere equilibrium Nu-2 = O(|Delta T|^(1/4)); mixed flow
        // similarly approaches its pure-forced coefficient (zero at rest).
        // The buoyant contribution to signed heat is o(Delta T), even when
        // dh/dTs diverges. Return the heat tangent, not a fictitious finite
        // coefficient derivative. Plate/cylinder cards retain their minima.
        if surface == bath
            && matches!(self.settings.law,
                SourceCellConvectionLaw::Natural {
                    geometry: NaturalConvectionGeometry::Sphere, ..
                } | SourceCellConvectionLaw::AidingVerticalPlate { .. })
        {
            if !energy.is_finite()
                || !(self.energy_bounds_j_kg[0]..=self.energy_bounds_j_kg[1]).contains(&energy)
                || !surface.is_finite()
                || !(self.settings.surface_temperature_bounds_k[0]
                    ..=self.settings.surface_temperature_bounds_k[1])
                    .contains(&surface)
            {
                return Err(DepositionError::MaterialDomain);
            }
            self.caloric.temperature_k(energy)?;
            let resolve = || -> Result<f64, LiquidExteriorConvectionError> {
                let sample = self.sample_transport(surface, false)?;
                let (beta, _) = self.expansion(surface, &sample, false)?;
                let length = self.settings.characteristic_length;
                let nu = match self.settings.law {
                    SourceCellConvectionLaw::Natural { gravity_m_s2, .. } => {
                        let ra = sample.state.rayleigh(gravity_m_s2, beta, 0.0, length.value())?;
                        NaturalConvectionGeometry::Sphere.evaluate(ra, sample.state.prandtl())?
                    }
                    SourceCellConvectionLaw::AidingVerticalPlate { settings } => {
                        let groups = settings.groups(sample.state, beta,
                            Temperature::new(surface), Temperature::new(bath), length)?;
                        match fs_convection::evaluate(
                            CorrelationId::ManningQureshiAidingVerticalPlate,
                            CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl),
                        ) {
                            Ok(nu) => nu,
                            // The original card validates all groups, including
                            // Pr, before its dedicated zero-driver refusal.
                            // Only at EXACT thermal equilibrium do heat and its
                            // tangent have the zero resting-flow limit. There
                            // is still no Nusselt evidence for this state.
                            Err(CorrelationError::InvalidGroup { axis: "max(Re,Ra)", .. })
                                if groups.reynolds == 0.0 && groups.rayleigh == 0.0 => return Ok(0.0),
                            Err(error) => return Err(error.into()),
                        }
                    }
                    _ => unreachable!("equilibrium branch only admits sphere and mixed plate"),
                };
                Ok(nu
                    .heat_transfer_coefficient(
                        ThermalConductivity::new(sample.state.thermal_conductivity_w_m_k()),
                        length,
                    )?
                    .value
                    .value())
            };
            let htc = resolve().map_err(|_| DepositionError::MaterialDomain)?;
            return Ok(ImplicitGrayExteriorFluxState {
                coefficients: ImplicitGraySurfaceState {
                    emissivity: self.settings.emissivity,
                    emissivity_enthalpy_derivative_kg_j: 0.0,
                    convection_w_m2_k: htc,
                    convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
                },
                emissivity_temperature_derivative_per_k: 0.0,
                convection_temperature_tangent_w_m2_k: htc,
            });
        }
        let state = self.state_at_surface(energy, surface, bath)?;
        let tangent = state.coefficients.convection_w_m2_k
            + (surface - bath) * state.convection_temperature_derivative_w_m2_k2;
        if !tangent.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitGrayExteriorFluxState {
            coefficients: state.coefficients,
            emissivity_temperature_derivative_per_k: state.emissivity_temperature_derivative_per_k,
            convection_temperature_tangent_w_m2_k: tangent,
        })
    }

    fn state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        if !energy.is_finite()
            || !(self.energy_bounds_j_kg[0]..=self.energy_bounds_j_kg[1]).contains(&energy)
        {
            return Err(DepositionError::MaterialDomain);
        }
        self.caloric.temperature_k(energy)?;
        let response = self
            .evaluate(surface, bath)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitGrayExteriorState {
            coefficients: ImplicitGraySurfaceState {
                emissivity: self.settings.emissivity,
                emissivity_enthalpy_derivative_kg_j: 0.0,
                convection_w_m2_k: response.coefficient.value.value(),
                convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
            },
            emissivity_temperature_derivative_per_k: 0.0,
            convection_temperature_derivative_w_m2_k2: response
                .coefficient_surface_derivative_w_m2_k2,
        })
    }
}

fn require_partial(
    evaluation: &NusseltEvaluation,
    axis: &'static str,
    varies: bool,
) -> Result<(), LiquidExteriorConvectionError> {
    if varies && !evaluation.evidence().sensitivity.d_qoi.contains_key(axis) {
        Err(LiquidExteriorConvectionError::MissingDerivative { axis })
    } else {
        Ok(())
    }
}
fn partial(
    evaluation: &NusseltEvaluation,
    axis: &'static str,
    input_derivative: f64,
) -> Result<f64, LiquidExteriorConvectionError> {
    if input_derivative == 0.0 {
        return Ok(0.0);
    }
    let derivative = evaluation
        .evidence()
        .sensitivity
        .d_qoi
        .get(axis)
        .ok_or(LiquidExteriorConvectionError::MissingDerivative { axis })?;
    finite(derivative * input_derivative)
}
fn finite(value: f64) -> Result<f64, LiquidExteriorConvectionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LiquidExteriorConvectionError::Invalid(
            "nonfinite source or tangent arithmetic",
        ))
    }
}
fn positive_interval([lo, hi]: [f64; 2]) -> bool {
    lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo
}
