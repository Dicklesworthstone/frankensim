//! Source-pinned property temperatures on the original finite gas/body heat residual.
//! Property evaluation has no independent energy or mass; bulk gas keeps both.
mod source;
pub use source::{FilmGasSource, PengRobinsonFilmSource};

use super::super::{
    AidingVerticalPlateSettings, evaluate_aiding_vertical_plate_from_material,
    evaluate_natural_from_material, evaluate_flat_plate_section_from_material,
    evaluate_turbulent_tube_section_from_material, mixed,
};
use super::*;
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// The existing correlation/geometry, with an explicit property-temperature convention.
#[derive(Clone, Copy, Debug)]
pub enum FilmGasHeatLaw {
    /// External flat-plate or cross-flow cylinder law at (T_wall + T_bulk)/2.
    /// Bulk-property duct laws must use `BulkForced`.
    Forced(ForcedGasSurface),
    /// Isothermal circular/rectangular duct law at the actual bulk temperature.
    /// This admits the existing laminar CWT, developing CWT, Dittus-Boelter and
    /// fully developed/developing Gnielinski cards. Their original geometry,
    /// regime and modest-property-variation assumptions still apply;
    /// no wall-viscosity correction is added.
    /// The source bank is checked even though the requested temperature is bulk.
    BulkForced(ForcedGasSurface),
    /// A physical streamwise strip; its boundary layer begins at the original
    /// upstream edge, not at this thermal patch's start.
    ForcedPlateSection(ForcedPlateGasSurface),
    /// Developing turbulent circular-tube section at the actual bulk state.
    /// The original physical entrance survives thermal-patch subdivision.
    /// This is the shared constant-property Gnielinski specialization K=1,
    /// with no wall-viscosity query or correction.
    BulkForcedTubeSection(ForcedTurbulentTubeSectionGasSurface),
    Natural(NaturalGasSurface),
    /// Laminar aiding mixed convection with explicitly bounded density departure.
    AidingVerticalPlate(AidingVerticalPlateGasSurface),
}

/// One flat-plate strip on the original finite gas/body heat residual.
///
/// Both distances are measured downstream from the physical leading edge.
/// Reversing the stream requires rebinding them from the opposite edge.
/// Geometry and area describe a real exposed patch; they are not inferred from
/// a cell centre. At common wall/film states adjacent strips recover the same
/// whole-plate heat. Local property refresh is a held-property approximation,
/// not a solution of a spatially non-isothermal boundary layer.
#[derive(Clone, Copy, Debug)]
pub struct ForcedPlateGasSurface {
    pub area_m2: f64,
    pub start: Length,
    pub end: Length,
    pub speed_m_s: f64,
    pub correlation: CorrelationId,
}

impl ForcedPlateGasSurface {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        let start = self.start.value();
        let end = self.end.value();
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !self.speed_m_s.is_finite() || self.speed_m_s < 0.0
            || (self.area_m2 > 0.0 && self.speed_m_s == 0.0)
            || !start.is_finite() || start < 0.0
            || !end.is_finite() || end <= start
            || (start > 0.0 && start / end == 0.0)
        {
            return Err(invalid(
                "finite gas plate sections need nonnegative area, positive speed for active heat transfer and finite 0 <= start < end from one leading edge",
            ));
        }
        if !matches!(self.correlation,
            CorrelationId::FlatPlateLaminarAverage | CorrelationId::FlatPlateTurbulentAverage)
        {
            return Err(invalid("finite gas plate sections require an original flat-plate mean correlation"));
        }
        Ok(())
    }
}

/// One inside-wall section of a smooth, straight turbulent circular tube.
///
/// `start` and `end` are distances from the same physical heated entrance,
/// not distances from this thermal patch. Flow reversal requires rebinding
/// them from the other entrance. The explicit area cannot exceed the actual
/// circular inside-wall section area pi D (end-start). Uniform bulk properties
/// and thermal driving make adjacent sections additive; local property refresh
/// does not solve spatially nonuniform entrance history. The original card's
/// Re, Pr and positive-endpoint L/D bounds remain mandatory on every trial.
#[derive(Clone, Copy, Debug)]
pub struct ForcedTurbulentTubeSectionGasSurface {
    pub area_m2: f64,
    pub diameter: Length,
    pub start: Length,
    pub end: Length,
    pub speed_m_s: f64,
}

impl ForcedTurbulentTubeSectionGasSurface {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        let diameter = self.diameter.value();
        let start = self.start.value();
        let end = self.end.value();
        let start_ratio = start / diameter;
        let end_ratio = end / diameter;
        let inside_area = core::f64::consts::PI * diameter * (end - start);
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !diameter.is_finite() || diameter <= 0.0
            || !start.is_finite() || start < 0.0
            || !end.is_finite() || end <= start
            || !start_ratio.is_finite() || (start > 0.0 && start_ratio == 0.0)
            || !end_ratio.is_finite() || end_ratio <= start_ratio
            || !inside_area.is_finite() || inside_area <= 0.0
            || self.area_m2 > inside_area
            || !self.speed_m_s.is_finite() || self.speed_m_s < 0.0
            || (self.area_m2 > 0.0 && self.speed_m_s == 0.0)
        {
            return Err(invalid(
                "finite gas turbulent tube sections need positive diameter, finite 0 <= start < end from one entrance, positive active speed and a patch no larger than pi D (end-start)",
            ));
        }
        Ok(())
    }
}

/// One isothermal vertical plate patch in a prescribed aiding free stream.
/// Its mean-film gas is a property query only; the bulk gas/body keep all storage.
/// The complete heat-solve bracket must remain in the aiding and density-valid
/// domain, including trial endpoints that the existing solver evaluates first.
#[derive(Clone, Copy, Debug)]
pub struct AidingVerticalPlateGasSurface {
    pub area_m2: f64,
    /// Plate height along the vertical free stream [m].
    pub characteristic_length: Length,
    pub settings: AidingVerticalPlateSettings,
}

impl AidingVerticalPlateGasSurface {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        self.settings.validate()?;
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
        {
            return Err(invalid("aiding mixed gas convection needs nonnegative area and positive plate height"));
        }
        Ok(())
    }
}

/// Source bank and constitutive boundary law for one finite gas/body patch.
/// The bank may contain more species than the current gas, or different initial
/// fractions. Each query reconstructs the ACTUAL active gas composition and
/// preserves its original source claims before evaluating the selected properties.
#[derive(Debug)]
pub struct FilmGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    pub model: &'a M,
    pub law: FilmGasHeatLaw,
}
impl<M: FilmGasSource> Copy for FilmGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for FilmGasSurface<'_, M> {
    fn clone(&self) -> Self { *self }
}

/// The actual bulk gas/body states and the separate property evaluation.
/// Bulk gas remains the thermodynamic/acoustic state; the property query never replaces it.
#[derive(Debug)]
pub struct FilmGasConvectionSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> {
    pub interface: FiniteGasConvectionSample<S, W>,
    /// Source-validated bulk properties for `BulkForced` and
    /// `BulkForcedTubeSection`, mean-film properties
    /// otherwise. None only for an explicitly isolated zero-area patch.
    pub property_gas: Option<S>,
}

impl FilmGasHeatLaw {
    pub(super) fn area(self) -> Result<f64, MaterialConvectionError> {
        match self {
            Self::Forced(surface) => {
                surface.validate()?;
                if surface.area_m2 > 0.0 && !matches!(surface.correlation,
                    CorrelationId::FlatPlateLaminarAverage
                    | CorrelationId::FlatPlateTurbulentAverage
                    | CorrelationId::ChurchillBernsteinCylinder)
                {
                    return Err(invalid(
                        "mean-film forced convection requires an external plate or cylinder correlation; use BulkForced for bulk-property duct laws",
                    ));
                }
                Ok(surface.area_m2)
            }
            Self::BulkForced(surface) => {
                surface.validate()?;
                if surface.area_m2 > 0.0 && !matches!(surface.correlation,
                    CorrelationId::CircularDuctLaminarCwt
                    | CorrelationId::CircularDuctHausen
                    | CorrelationId::RectangularDuctLaminarCwt
                    | CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
                    | CorrelationId::DittusBoelter
                    | CorrelationId::Gnielinski
                    | CorrelationId::GnielinskiDeveloping)
                {
                    return Err(invalid(
                        "bulk forced convection requires an isothermal duct correlation without a separate wall-viscosity input",
                    ));
                }
                Ok(surface.area_m2)
            }
            Self::BulkForcedTubeSection(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
            Self::ForcedPlateSection(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
            Self::Natural(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
            Self::AidingVerticalPlate(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
        }
    }
}

impl<M: FilmGasSource> FilmGasSurface<'_, M> {
    /// Re-query source-pinned properties at the actual bulk pressure and the
    /// selected reference temperature: T_bulk for `BulkForced` and
    /// `BulkForcedTubeSection`, otherwise
    /// (T_wall + T_bulk)/2. Both signed heat flow and Rayleigh forcing retain
    /// the full wall-to-bulk temperature contrast. The property query is not
    /// a replacement gas endpoint or a second finite-capacity reservoir.
    ///
    /// The original correlation must support this property convention and the
    /// declared geometry. The original ideal or explicit-branch real-gas owner
    /// retains frozen composition and its single-phase source domain. This is
    /// a property query with no boundary-layer or chemical heat solver.
    /// The callback works with existing powered/reactive/multi-boundary owners;
    /// their complete exposed-area checks and source budgets remain mandatory.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FilmGasConvectionSample<M::State>>, MaterialConvectionError> {
        self.resolve_at_temperature(gas, body.temperature_k(), body)
    }

    fn resolve_at_temperature<W>(
        self,
        gas: &M::State,
        wall: f64,
        body: W,
    ) -> Result<ConductanceSample<FilmGasConvectionSample<M::State, W>>, MaterialConvectionError> {
        let area_m2 = self.law.area()?;
        let bulk = gas.temperature_k();
        if !wall.is_finite() || wall <= 0.0 {
            return Err(invalid("film convection needs a positive body temperature"));
        }
        if area_m2 == 0.0 {
            return Ok(ConductanceSample {
                conductance_w_k: 0.0,
                source: FilmGasConvectionSample {
                    interface: FiniteGasConvectionSample {
                        gas: gas.clone(),
                        body,
                        area_m2,
                        evaluation: None,
                    },
                    property_gas: None,
                },
            });
        }
        let property_temperature = match self.law {
            FilmGasHeatLaw::BulkForced(_) | FilmGasHeatLaw::BulkForcedTubeSection(_) => bulk,
            _ => 0.5 * wall + 0.5 * bulk,
        };
        // Re-query even at T_bulk: cloning the snapshot here would bypass the
        // original source-bank and real-gas branch-identity checks.
        let property = self
            .model
            .evaluate_property_temperature(gas, Temperature::new(property_temperature))?;
        let expansion = || -> Result<f64, MaterialConvectionError> {
            let thermal = property.thermal_transport_state()?;
            if thermal.transport() != property.transport_state()? || thermal.temperature_k() != property.temperature_k() {
                return Err(invalid("film thermal and momentum views must describe the same material state"));
            }
            thermal.volumetric_expansion_per_k()
                .ok_or_else(|| invalid("film buoyancy requires source-owned volumetric expansion"))
        };
        let evaluation = match self.law {
            FilmGasHeatLaw::ForcedPlateSection(surface) => {
                let evaluated = evaluate_flat_plate_section_from_material(
                    surface.correlation, &property, surface.speed_m_s, surface.start, surface.end,
                ).map_err(MaterialConvectionError::from)?;
                let (transport, nusselt, coefficient, section) = evaluated.into_parts();
                Some(GasConvectionEvaluation {
                    transport, nusselt, coefficient, section: Some(section), tube_section: None,
                })
            }
            FilmGasHeatLaw::BulkForcedTubeSection(surface) => {
                let evaluated = evaluate_turbulent_tube_section_from_material(
                    &property, surface.speed_m_s, surface.diameter, surface.start, surface.end,
                )?;
                let (transport, nusselt, coefficient, tube_section) = evaluated.into_parts();
                Some(GasConvectionEvaluation {
                    transport, nusselt, coefficient, section: None, tube_section: Some(tube_section),
                })
            }
            FilmGasHeatLaw::Forced(surface) | FilmGasHeatLaw::BulkForced(surface) => {
                let mut settings = surface.settings;
                if wall > bulk {
                    settings.direction = ThermalDirection::HeatingFluid;
                } else if wall < bulk {
                    settings.direction = ThermalDirection::CoolingFluid;
                }
                Some(GasConvectionEvaluation::from_material(
                    evaluate_forced_from_material(
                        surface.correlation,
                        &property,
                        surface.speed_m_s,
                        surface.characteristic_length,
                        settings,
                    )?,
                ))
            }
            FilmGasHeatLaw::Natural(surface) => {
                if wall == bulk {
                    None
                } else {
                    Some(GasConvectionEvaluation::from_material(
                        evaluate_natural_from_material(
                            surface.geometry,
                            &property,
                            surface.gravity_m_s2,
                            expansion()?,
                            Temperature::new(wall),
                            Temperature::new(bulk),
                            surface.characteristic_length,
                        )?,
                    ))
                }
            }
            FilmGasHeatLaw::AidingVerticalPlate(surface) => {
                if wall == bulk && surface.settings.speed_up_m_s == 0.0 {
                    mixed::inputs(property.transport_state()?, surface.settings,
                        expansion()?, Temperature::new(wall),
                        Temperature::new(bulk), surface.characteristic_length)?;
                    None
                } else {
                    Some(GasConvectionEvaluation::from_material(
                        evaluate_aiding_vertical_plate_from_material(
                            &property, surface.settings, expansion()?,
                            Temperature::new(wall), Temperature::new(bulk),
                            surface.characteristic_length,
                        )?,
                    ))
                }
            }
        };
        let conductance_w_k = match &evaluation {
            None => 0.0,
            Some(value) => {
                let conductance = area_m2 * value.coefficient.value.value();
                if !conductance.is_finite() || conductance <= 0.0 {
                    return Err(invalid(
                        "film convective conductance overflowed or vanished",
                    ));
                }
                conductance
            }
        };
        Ok(ConductanceSample {
            conductance_w_k,
            source: FilmGasConvectionSample {
                interface: FiniteGasConvectionSample {
                    gas: gas.clone(),
                    body,
                    area_m2,
                    evaluation,
                },
                property_gas: Some(property),
            },
        })
    }

    /// Advance heat, source-backed gas energy, finite body enthalpy and piston
    /// momentum together through the existing heat/motion residual. The same
    /// proposed bulk gas determines pressure work, storage and acoustic state;
    /// only the transfer coefficient uses the selected source-pinned properties.
    /// Both original physical endpoints are returned together or neither.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate<G: PistonGas<State = M::State>>(
        self,
        piston: &GasPiston<G>,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<FilmGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        let area = self.law.area().map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || area > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid(
                "film convective patch exceeds the finite body's exposed area",
            ));
        }
        piston.thermal_candidate(
            input,
            body,
            controls,
            exchange,
            |gas, wall| self.resolve(gas, wall),
            cancelled,
        )
    }

    /// Resolve source properties at the coupled wall and fixed-volume gas
    /// endpoints through the original reciprocal heat and wall-feedback solve.
    /// The retained wall query is a temperature, not an inferred phase fraction;
    /// the body owner retains its actual enthalpy, latent heat and phase state.
    ///
    /// Positive area must cover the complete exposed body surface so its local
    /// coefficient reaches the original Biot check without area averaging.
    /// Radiation remains external. All gas/heat and wall-temperature brackets
    /// must stay within the source and correlation domains, including mixed
    /// convection's aiding-flow and density-departure limits.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = M::State>>(
        self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<WallConvectionSample<FilmGasConvectionSample<M::State, Temperature>>, G>,
        >,
        ConvectionError<WallConvectionError<MaterialConvectionError>>,
    > {
        cx.checkpoint().map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        let area = self.law.area().map_err(|error| {
            ConvectionError::Convection(WallConvectionError::Source(error))
        })?;
        if !body.surface_area_m2.is_finite() || body.surface_area_m2 <= 0.0
            || (area != 0.0 && area != body.surface_area_m2)
        {
            return Err(ConvectionError::Thermal(ThermalPistonError::Invalid(
                "isochoric film convection must cover the complete body surface or be explicitly insulated",
            )));
        }
        environment.advance_with_wall_convection(cx, body, wall, |gas, temperature_k| {
            let sample = self.resolve_at_temperature(gas, temperature_k, Temperature::new(temperature_k))?;
            let coefficient_w_per_m2_k = sample.source.interface.evaluation.as_ref()
                .map_or(0.0, |value| value.coefficient.value.value());
            Ok(ConvectionCoefficient { coefficient_w_per_m2_k, source: sample.source })
        })
    }
}
