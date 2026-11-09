//! Prescribed stream mass flux on the original finite gas/body heat residual.
use super::*;

// Keep the original film-facing public path while all finite owners share
// one positive-input and original-EOS velocity implementation.
pub use super::super::GasMassFlowState;
pub(super) use super::super::mass_flow::MassFlowInput;

impl MassFlowInput {
    fn admit_geometry(self, law: FilmGasHeatLaw) -> Result<(), MaterialConvectionError> {
        let rectangular_area = |diameter: f64, alpha: f64| {
            if !(alpha > 0.0 && alpha <= 1.0) {
                return f64::NAN;
            }
            diameter * diameter * (1.0 + alpha) * (1.0 + alpha) / (4.0 * alpha)
        };
        let circular_area = |diameter: f64| core::f64::consts::PI * diameter * diameter / 4.0;
        let expected = match law {
            FilmGasHeatLaw::BulkForced(surface) => {
                let diameter = surface.characteristic_length.value();
                if matches!(
                    surface.correlation,
                    CorrelationId::RectangularDuctLaminarCwt
                        | CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
                ) {
                    rectangular_area(diameter, surface.settings.aspect_ratio.unwrap_or(f64::NAN))
                } else {
                    circular_area(diameter)
                }
            }
            FilmGasHeatLaw::BulkForcedTubeSection(surface)
            | FilmGasHeatLaw::BulkForcedLaminarTubeSection(surface) => {
                circular_area(surface.diameter.value())
            }
            FilmGasHeatLaw::BulkForcedRectangularDuctSection(surface) => {
                rectangular_area(surface.hydraulic_diameter.value(), surface.aspect_ratio)
            }
            _ => return Ok(()), // An external stream's cross section is authored independently.
        };
        self.admit_flow_area(expected)
    }

    pub(super) fn bind(
        self,
        law: FilmGasHeatLaw,
        gas: &impl PistonGasState,
    ) -> Result<(FilmGasHeatLaw, GasMassFlowState), MaterialConvectionError> {
        let state = self.resolve(gas)?;
        Ok((with_speed(law, state.speed_m_s)?, state))
    }
}

/// Prescribed stream/recirculation mass rate for an original forced heat law.
///
/// Each nonlinear trial uses `U = mdot/(rho_bulk A_flow)` with the original
/// bulk gas state, even when the correlation samples mean-film properties.
/// Pressure, composition, source claims and physical inlet coordinates remain
/// unchanged. Internal duct flow area must match the actual circular/rectangular
/// cross section; an external stream's flow area is independently authored.
/// A positive scalar rate has no flow direction; callers must bind
/// the physical entrance explicitly after any reversal. Natural and oriented
/// mixed convection are unsupported by this forced-speed selection.
///
/// This rate controls only the convection velocity. It does not add or remove
/// gas mass or enthalpy, prescribe fan/pump work, resolve axial advection, or
/// change either finite energy owner. Both native coupled solvers are shared
/// with `FilmGasSurface`; zero thermal area remains exact insulation.
#[derive(Debug)]
pub struct PrescribedMassFlowGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    surface: FilmGasSurface<'a, M>,
    flow: MassFlowInput,
}
impl<M: FilmGasSource> Copy for PrescribedMassFlowGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for PrescribedMassFlowGasSurface<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, M: FilmGasSource> FilmGasSurface<'a, M> {
    /// Replace the authored fixed speed with a positive prescribed stream mass
    /// rate and cross section, using each actual trial's original bulk density.
    pub fn with_prescribed_mass_flow(
        self,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<PrescribedMassFlowGasSurface<'a, M>, MaterialConvectionError> {
        let flow = MassFlowInput::new(mass_flow_kg_s, flow_area_m2)?;
        // Only classify the law here; no reference state is invented or queried.
        with_speed(self.law, 0.0)?;
        flow.admit_geometry(self.law)?;
        Ok(PrescribedMassFlowGasSurface {
            surface: self,
            flow,
        })
    }
}

impl<M: FilmGasSource> PrescribedMassFlowGasSurface<'_, M> {
    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.flow.mass_flow_kg_s()
    }
    pub const fn flow_area_m2(&self) -> f64 {
        self.flow.flow_area_m2()
    }

    /// Resolve the original source and current bulk-density velocity together.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FilmGasConvectionSample<M::State>>, MaterialConvectionError> {
        self.surface
            .resolve_at_temperature(Some(self.flow), gas, body.temperature_k(), body)
    }

    /// Resolve a single- or shared-body wall trial with the original bulk-density
    /// mass-flow conversion. Film density never determines stream velocity.
    /// Complete exposed-area admission and the coefficient/source meaning are
    /// shared with [`FilmGasSurface::resolve_isochoric_coefficient`].
    pub fn resolve_isochoric_coefficient(
        self,
        gas: &M::State,
        body_surface_area_m2: f64,
        wall_temperature: Temperature,
    ) -> Result<
        ConvectionCoefficient<FilmGasConvectionSample<M::State, Temperature>>,
        MaterialConvectionError,
    > {
        self.surface.resolve_isochoric_coefficient_with(
            Some(self.flow), gas, body_surface_area_m2, wall_temperature,
        )
    }

    /// Propose the same reciprocal gas/body heat and piston motion step.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate<G: PistonGas<State = M::State>>(
        self,
        piston: &GasPiston<G>,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<FilmGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        self.surface.finite_body_candidate_with(
            Some(self.flow),
            piston,
            input,
            body,
            controls,
            exchange,
            cancelled,
        )
    }

    /// Propose the same fixed-volume gas/body and wall-feedback step.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = M::State>>(
        self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<
                WallConvectionSample<FilmGasConvectionSample<M::State, Temperature>>,
                G,
            >,
        >,
        ConvectionError<WallConvectionError<MaterialConvectionError>>,
    > {
        self.surface
            .isochoric_body_candidate_with(Some(self.flow), cx, environment, body, wall)
    }
}

impl<'g, M: FilmGasSource> PrescribedMassFlowGasSurface<'g, M> {
    /// Retain original bulk-density mass flux while resolving the conducting
    /// wall and source film at their actual interface temperature.
    pub fn with_conducting_layer<'a, 'm>(
        self,
        layer: &'a fs_material::thermal::conduction::source::ResolvedConductionLayer<'m>,
        geometry: fs_material::thermal::conduction::ReservoirLayerGeometry,
        config: fs_material::thermal::conduction::LayerConductionConfig,
    ) -> Result<
        super::super::conjugate::ConjugateFilmGasSurface<'g, 'a, 'm, M>,
        super::super::conjugate::ConjugateGasError,
    > {
        super::super::conjugate::ConjugateFilmGasSurface::from_film(
            self.surface, Some(self.flow), layer, geometry, config,
        )
    }
}

fn with_speed(
    law: FilmGasHeatLaw,
    speed_m_s: f64,
) -> Result<FilmGasHeatLaw, MaterialConvectionError> {
    Ok(match law {
        FilmGasHeatLaw::Forced(surface) => FilmGasHeatLaw::Forced(ForcedGasSurface {
            speed_m_s,
            ..surface
        }),
        FilmGasHeatLaw::BulkForced(surface) => FilmGasHeatLaw::BulkForced(ForcedGasSurface {
            speed_m_s,
            ..surface
        }),
        FilmGasHeatLaw::ForcedPlateSection(surface) => {
            FilmGasHeatLaw::ForcedPlateSection(ForcedPlateGasSurface {
                speed_m_s,
                ..surface
            })
        }
        FilmGasHeatLaw::BulkForcedTubeSection(surface) => {
            FilmGasHeatLaw::BulkForcedTubeSection(ForcedTubeSectionGasSurface {
                speed_m_s,
                ..surface
            })
        }
        FilmGasHeatLaw::BulkForcedLaminarTubeSection(surface) => {
            FilmGasHeatLaw::BulkForcedLaminarTubeSection(ForcedTubeSectionGasSurface {
                speed_m_s,
                ..surface
            })
        }
        FilmGasHeatLaw::BulkForcedRectangularDuctSection(surface) => {
            FilmGasHeatLaw::BulkForcedRectangularDuctSection(
                ForcedRectangularDuctSectionGasSurface {
                    speed_m_s,
                    ..surface
                },
            )
        }
        FilmGasHeatLaw::Natural(_) | FilmGasHeatLaw::AidingVerticalPlate(_) => {
            return Err(invalid(
                "prescribed mass flow requires a forced heat law; natural or oriented mixed flow is not inferred",
            ));
        }
    })
}
