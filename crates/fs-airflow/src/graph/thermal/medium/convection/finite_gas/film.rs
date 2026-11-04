//! Film-temperature properties on the original finite gas/body heat residual.
//! Film evaluation has no independent energy or mass; bulk gas keeps both.
mod source;
pub use source::{FilmGasSource, PengRobinsonFilmSource};

use super::super::{
    AidingVerticalPlateSettings, evaluate_aiding_vertical_plate_from_material,
    evaluate_natural_from_material, evaluate_flat_plate_section_from_material, mixed,
};
use super::*;
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// The existing correlation/geometry, with an explicit mean-film convention.
#[derive(Clone, Copy, Debug)]
pub enum FilmGasHeatLaw {
    Forced(ForcedGasSurface),
    /// A physical streamwise strip; its boundary layer begins at the original
    /// upstream edge, not at this thermal patch's start.
    ForcedPlateSection(ForcedPlateGasSurface),
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
            || !start.is_finite() || start < 0.0
            || !end.is_finite() || end <= start
            || (start > 0.0 && start / end == 0.0)
        {
            return Err(invalid(
                "finite gas plate sections need nonnegative area/speed and finite 0 <= start < end from one leading edge",
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
/// preserves its original source claims before evaluating film properties.
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
/// Bulk gas remains the thermodynamic/acoustic state; film never replaces it.
#[derive(Debug)]
pub struct FilmGasConvectionSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> {
    pub interface: FiniteGasConvectionSample<S, W>,
    /// None only for an explicitly isolated zero-area patch.
    pub property_gas: Option<S>,
}

impl FilmGasHeatLaw {
    pub(super) fn area(self) -> Result<f64, MaterialConvectionError> {
        match self {
            Self::Forced(surface) => {
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
    /// Re-query source-pinned properties at (T_wall + T_bulk)/2 and actual bulk
    /// pressure. Both signed heat flow and Rayleigh forcing retain the FULL
    /// wall-to-bulk temperature contrast. The film is not a replacement gas
    /// endpoint, a second parcel, or an invented finite-capacity reservoir.
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
        let film = self
            .model
            .evaluate_property_temperature(gas, Temperature::new(0.5 * wall + 0.5 * bulk))?;
        let expansion = || -> Result<f64, MaterialConvectionError> {
            let thermal = film.thermal_transport_state()?;
            if thermal.transport() != film.transport_state()? || thermal.temperature_k() != film.temperature_k() {
                return Err(invalid("film thermal and momentum views must describe the same material state"));
            }
            thermal.volumetric_expansion_per_k()
                .ok_or_else(|| invalid("film buoyancy requires source-owned volumetric expansion"))
        };
        let evaluation = match self.law {
            FilmGasHeatLaw::ForcedPlateSection(surface) => {
                let evaluated = evaluate_flat_plate_section_from_material(
                    surface.correlation, &film, surface.speed_m_s, surface.start, surface.end,
                ).map_err(MaterialConvectionError::from)?;
                let (transport, nusselt, coefficient, section) = evaluated.into_parts();
                Some(GasConvectionEvaluation {
                    transport, nusselt, coefficient, section: Some(section),
                })
            }
            FilmGasHeatLaw::Forced(surface) => {
                let mut settings = surface.settings;
                if wall > bulk {
                    settings.direction = ThermalDirection::HeatingFluid;
                } else if wall < bulk {
                    settings.direction = ThermalDirection::CoolingFluid;
                }
                Some(GasConvectionEvaluation::from_material(
                    evaluate_forced_from_material(
                        surface.correlation,
                        &film,
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
                            &film,
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
                    mixed::inputs(film.transport_state()?, surface.settings,
                        expansion()?, Temperature::new(wall),
                        Temperature::new(bulk), surface.characteristic_length)?;
                    None
                } else {
                    Some(GasConvectionEvaluation::from_material(
                        evaluate_aiding_vertical_plate_from_material(
                            &film, surface.settings, expansion()?,
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
                property_gas: Some(film),
            },
        })
    }

    /// Advance heat, source-backed gas energy, finite body enthalpy and piston
    /// momentum together through the existing heat/motion residual. The same
    /// proposed bulk gas determines pressure work, storage and acoustic state;
    /// only the transfer coefficient uses its source-pinned film properties.
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

    /// Resolve film properties at the coupled wall and fixed-volume gas
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
        if area > 0.0 && matches!(self.law, FilmGasHeatLaw::Forced(ForcedGasSurface {
            correlation: CorrelationId::CircularDuctLaminarChf | CorrelationId::RectangularDuctLaminarChf,
            ..
        })) {
            return Err(ConvectionError::Convection(WallConvectionError::Source(invalid(
                "an isothermal finite body cannot impose a constant-heat-flux duct correlation",
            ))));
        }
        environment.advance_with_wall_convection(cx, body, wall, |gas, temperature_k| {
            let sample = self.resolve_at_temperature(gas, temperature_k, Temperature::new(temperature_k))?;
            let coefficient_w_per_m2_k = sample.source.interface.evaluation.as_ref()
                .map_or(0.0, |value| value.coefficient.value.value());
            Ok(ConvectionCoefficient { coefficient_w_per_m2_k, source: sample.source })
        })
    }
}
