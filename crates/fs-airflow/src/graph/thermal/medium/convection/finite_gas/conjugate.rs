//! Original conducting-wall surfaces on the existing finite gas/body owners.
//!
//! The material layer eliminates its actual fluid-facing temperature against
//! the original gas film. Only the combined conductance enters the reciprocal
//! heat solve; film hA is diagnostic and is never applied as a second source.

use super::film::{FilmGasConvectionSample, FilmGasSource, FilmGasSurface};
use super::mass_flow::MassFlowInput;
use super::*;
use fs_matdb::PropertyLineSlice;
use fs_material::thermal::conduction::{
    ConvectiveLayerError, ConvectiveLayerResponse, LayerConductionConfig, ReservoirLayerGeometry,
    source::{LayerSourceError, ResolvedConductionLayer},
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;
use fs_thermal::enthalpy::DepositionError;
use std::{cell::RefCell, sync::Arc};

/// Original film/material failures or an unclosed physical surface balance.
#[derive(Debug)]
pub enum ConjugateGasError {
    Invalid(&'static str),
    Film(MaterialConvectionError),
    Layer(LayerSourceError),
    Thermal(DepositionError),
    SurfaceBalance { residual_w: f64, iterations: usize },
}
impl core::fmt::Display for ConjugateGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas conducting surface: {self:?}")
    }
}
impl std::error::Error for ConjugateGasError {}
impl From<MaterialConvectionError> for ConjugateGasError {
    fn from(error: MaterialConvectionError) -> Self {
        Self::Film(error)
    }
}
impl From<LayerSourceError> for ConjugateGasError {
    fn from(error: LayerSourceError) -> Self {
        Self::Layer(error)
    }
}
impl From<DepositionError> for ConjugateGasError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// The original body state and the film at its distinct, eliminated surface.
/// The layer and insulated_source fields are mutually exclusive. Zero film
/// area retains its original sample without evaluating a conducting root.
#[derive(Debug)]
pub struct ConjugateGasConvectionSample<F, W = EquilibriumPhaseState> {
    pub body: W,
    pub body_temperature_k: f64,
    pub bulk_temperature_k: f64,
    pub area_m2: f64,
    pub geometry: ReservoirLayerGeometry,
    pub conductivity_source: Arc<PropertyLineSlice>,
    pub layer: Option<ConvectiveLayerResponse<F>>,
    pub insulated_source: Option<F>,
    /// Combined G times the BODY-to-bulk difference, applied only once.
    pub applied_heat_into_gas_w: f64,
    /// Applied heat minus the original material owner's accepted film power.
    pub representation_error_w: f64,
}

pub type ConjugateFilmGasSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> =
    ConjugateGasConvectionSample<FilmGasConvectionSample<S, Temperature>, W>;

pub(super) mod sealed {
    pub trait Film {}
}

/// Original source film admitted at its actual fluid-facing temperature.
///
/// Implemented only by the existing source-pinned film, Sieder--Tate tube and
/// Whitaker sphere owners. Their exact property convention, geometry, selected
/// regime and prescribed-flow inputs remain authoritative. This seam supplies
/// a coefficient and original source payload, never storage or another heat solve.
pub trait ConjugateGasFilm: sealed::Film + Copy {
    type State: PistonGasState + ThermalFluidTransport;
    type Sample;

    /// Admit original geometry and flow; zero explicitly insulates the patch.
    fn admitted_area_m2(&self, gas: &Self::State) -> Result<f64, MaterialConvectionError>;

    /// Resolve only the original film at the material owner's actual surface.
    fn resolve_film_at_surface(
        &self,
        gas: &Self::State,
        surface: Temperature,
    ) -> Result<ConductanceSample<Self::Sample>, MaterialConvectionError>;

    /// Add the original source-integrated conducting path and shared native
    /// finite-body transactions without replacing this film or its flow law.
    fn with_conducting_layer<'a, 'm>(
        self,
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ReservoirLayerGeometry,
        config: LayerConductionConfig,
    ) -> Result<ConjugateGasSurface<'a, 'm, Self>, ConjugateGasError> {
        ConjugateGasSurface::new(self, layer, geometry, config)
    }
}

/// Opaque fixed-speed or prescribed-flow selection behind the original alias.
/// Its fields stay private; construction remains with the original film API.
#[doc(hidden)]
#[derive(Debug)]
pub struct ConjugateFilmSelection<'g, M: FilmGasSource = ThermochemicalGasModel> {
    film: FilmGasSurface<'g, M>,
    flow: Option<MassFlowInput>,
}
impl<M: FilmGasSource> Copy for ConjugateFilmSelection<'_, M> {}
impl<M: FilmGasSource> Clone for ConjugateFilmSelection<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<M: FilmGasSource> sealed::Film for ConjugateFilmSelection<'_, M> {}
impl<M: FilmGasSource> ConjugateGasFilm for ConjugateFilmSelection<'_, M> {
    type State = M::State;
    type Sample = FilmGasConvectionSample<M::State, Temperature>;

    fn admitted_area_m2(&self, gas: &Self::State) -> Result<f64, MaterialConvectionError> {
        let (law, _) = self.film.bind_flow(self.flow, gas)?;
        law.area()
    }

    fn resolve_film_at_surface(
        &self,
        gas: &Self::State,
        surface: Temperature,
    ) -> Result<ConductanceSample<Self::Sample>, MaterialConvectionError> {
        self.film
            .resolve_at_temperature(self.flow, gas, surface.value(), surface)
    }
}

/// Original type name shared by BOTH fixed-speed and prescribed-film constructors.
/// The same typed film payload is retained by the generic conducting owner.
pub type ConjugateFilmGasSurface<'g, 'a, 'm, M = ThermochemicalGasModel> =
    ConjugateGasSurface<'a, 'm, ConjugateFilmSelection<'g, M>>;

#[derive(Clone, Copy, Debug)]
struct ConductingPath<'a, 'm> {
    layer: &'a ResolvedConductionLayer<'m>,
    geometry: ReservoirLayerGeometry,
    config: LayerConductionConfig,
}

impl ConductingPath<'_, '_> {
    fn admit_state(&self, state: EquilibriumPhaseState) -> Result<(), ConjugateGasError> {
        if state.phase_curve_identity() != self.layer.caloric_curve().identity() {
            return Err(ConjugateGasError::Invalid(
                "conducting surface requires the original body's caloric identity",
            ));
        }
        let supported = self
            .layer
            .temperature_at_enthalpy(state.specific_enthalpy_j_kg())?;
        // The enthalpy gate may snap a source endpoint after caloric inversion.
        // Final piston states retain their actual native T without a film query.
        if state.temperature_k() != supported {
            self.layer
                .conductivity_curve()
                .conductivity_w_m_k(state.temperature_k())?;
        }
        Ok(())
    }

    fn admit_input(&self, body: UniformEnthalpyStepInput<'_>) -> Result<(), ConjugateGasError> {
        if !core::ptr::eq(self.layer.caloric_curve(), body.curve) {
            return Err(ConjugateGasError::Invalid(
                "conducting path and finite body must share the exact caloric object",
            ));
        }
        self.admit_state(body.initial)
    }

    // The original source callback owns its law and concrete payload. Every
    // supported film reuses this SAME material solve and native heat transaction.
    #[allow(clippy::too_many_arguments)]
    fn resolve_with<F, W>(
        &self,
        body_temperature_k: f64,
        bulk_temperature_k: f64,
        area_m2: f64,
        body: W,
        mut film: impl FnMut(Temperature) -> Result<ConductanceSample<F>, MaterialConvectionError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<ConjugateGasConvectionSample<F, W>>, ConjugateGasError> {
        checkpoint()?;
        if !body_temperature_k.is_finite()
            || body_temperature_k <= 0.0
            || !bulk_temperature_k.is_finite()
            || bulk_temperature_k <= 0.0
        {
            return Err(ConjugateGasError::Invalid(
                "positive finite body and bulk temperatures required",
            ));
        }
        let (conductance_w_k, layer, insulated_source, applied, representation_error_w) =
            if area_m2 == 0.0 {
                let zero = film(Temperature::new(body_temperature_k))?;
                if zero.conductance_w_k != 0.0 {
                    return Err(ConjugateGasError::Invalid(
                        "zero-area film must remain exactly insulated",
                    ));
                }
                (0.0, None, Some(zero.source), 0.0, 0.0)
            } else {
                let layer = self
                    .geometry
                    .evaluate_convective_material_to_bath(
                        self.layer.conductivity_curve(),
                        [body_temperature_k, bulk_temperature_k],
                        self.config,
                        |surface| {
                            let sample = film(Temperature::new(surface))?;
                            let power = sample.conductance_w_k * (surface - bulk_temperature_k);
                            Ok((power, sample.conductance_w_k, sample.source))
                        },
                        &mut checkpoint,
                    )
                    .map_err(|error| match error {
                        ConvectiveLayerError::Thermal(error) => ConjugateGasError::Thermal(error),
                        ConvectiveLayerError::Film(error) => ConjugateGasError::Film(error),
                        ConvectiveLayerError::SurfaceBalance {
                            residual_w,
                            iterations,
                        } => ConjugateGasError::SurfaceBalance {
                            residual_w,
                            iterations,
                        },
                    })?;
                let difference = body_temperature_k - bulk_temperature_k;
                let applied = layer.conductance_w_k * difference;
                let error = applied - layer.power_to_bath_w;
                let allowance = 8.0 * f64::EPSILON * applied.abs().max(layer.power_to_bath_w.abs());
                if !applied.is_finite()
                    || !error.is_finite()
                    || error.abs() > allowance
                    || (difference != 0.0 && applied == 0.0)
                    || (difference == 0.0 && layer.power_to_bath_w != 0.0)
                {
                    return Err(DepositionError::Arithmetic.into());
                }
                (layer.conductance_w_k, Some(layer), None, applied, error)
            };
        checkpoint()?;
        Ok(ConductanceSample {
            conductance_w_k,
            source: ConjugateGasConvectionSample {
                body,
                body_temperature_k,
                bulk_temperature_k,
                area_m2,
                geometry: self.geometry,
                conductivity_source: Arc::clone(self.layer.source_line()),
                layer,
                insulated_source,
                applied_heat_into_gas_w: applied,
                representation_error_w,
            },
        })
    }
}

impl<F, W> ConjugateGasConvectionSample<F, W> {
    fn admit_work(&self, duration_s: f64, tolerance_j: f64) -> Result<(), ConjugateGasError> {
        if let Some(layer) = &self.layer {
            let residual_j = duration_s * layer.heat_residual_w;
            if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
                return Err(ConjugateGasError::SurfaceBalance {
                    residual_w: layer.heat_residual_w,
                    iterations: layer.iterations,
                });
            }
        }
        Ok(())
    }
}

/// Source-integrated pure-phase path from the stored body trace to its film.
///
/// Construct from an original source film's with_conducting_layer selection.
/// The path uses the SAME original body caloric curve and supported
/// conductivity source. Its physical area must equal the positive film area.
/// Zero film area retains prepared geometry/caloric bindings and participating
/// body-state admission, but skips conductivity integration and film properties.
///
/// This is a homogeneous, source-free, quasi-steady normal path with no extra
/// thermal storage or spreading model. A uniform body trace and well-mixed gas
/// remain reduced-model assumptions; this does not resolve transient wall fields.
/// Film geometry, physical entrance coordinates, source and regime limits apply
/// at every actual surface. Prescribed flow still uses the original BULK density.
///
/// Isochoric use requires disabled radiation: the original radiation owner acts
/// at body temperature and cannot silently be moved to this distinct surface.
/// Its native Biot check still acts on combined G/area and its original body
/// conductivity; the material path is not a waiver of uniform-core admission.
#[derive(Debug)]
pub struct ConjugateGasSurface<'a, 'm, F: ConjugateGasFilm> {
    film: F,
    path: ConductingPath<'a, 'm>,
}
impl<F: ConjugateGasFilm> Copy for ConjugateGasSurface<'_, '_, F> {}
impl<F: ConjugateGasFilm> Clone for ConjugateGasSurface<'_, '_, F> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, 'm, F: ConjugateGasFilm> ConjugateGasSurface<'a, 'm, F> {
    /// Bind an original admitted source film to the same material path.
    pub fn new(
        film: F,
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ReservoirLayerGeometry,
        config: LayerConductionConfig,
    ) -> Result<Self, ConjugateGasError> {
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(ConjugateGasError::Invalid(
                "finite nonnegative layer residual budgets and positive work bound required",
            ));
        }
        Ok(Self {
            film,
            path: ConductingPath {
                layer,
                geometry,
                config,
            },
        })
    }

    fn area(&self, gas: &F::State) -> Result<f64, ConjugateGasError> {
        let area = self.film.admitted_area_m2(gas)?;
        if area != 0.0 && area != self.path.geometry.surface_area_m2() {
            return Err(ConjugateGasError::Invalid(
                "positive film and physical conducting face areas must match",
            ));
        }
        Ok(area)
    }

    fn resolve_at_temperature<W>(
        &self,
        gas: &F::State,
        temperature_k: f64,
        body: W,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<ConjugateGasConvectionSample<F::Sample, W>>, ConjugateGasError>
    {
        let area = self.area(gas)?;
        self.path.resolve_with(
            temperature_k,
            gas.temperature_k(),
            area,
            body,
            |surface| self.film.resolve_film_at_surface(gas, surface),
            checkpoint,
        )
    }

    /// Resolve actual surface Ts and combined conductance without advancing heat.
    /// The actual body state must remain on the prepared original pure phase.
    pub fn resolve(
        &self,
        gas: &F::State,
        body: EquilibriumPhaseState,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<ConjugateGasConvectionSample<F::Sample>>, ConjugateGasError> {
        self.path.admit_state(body)?;
        self.resolve_at_temperature(gas, body.temperature_k(), body, checkpoint)
    }

    /// Use the original midpoint gas/body/piston transaction and one reciprocal Q.
    /// Layer residual watts must also fit the enclosing source budget in joules.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate<G: PistonGas<State = F::State>>(
        &self,
        piston: &GasPiston<G>,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<ConjugateGasConvectionSample<F::Sample>, G>,
        ThermalPistonError<ConjugateGasError>,
    > {
        let cancelled = RefCell::new(cancelled);
        if (cancelled.borrow_mut())() {
            return Err(ThermalPistonError::Cancelled);
        }
        self.path
            .admit_input(body)
            .map_err(ThermalPistonError::Interface)?;
        let area = self
            .area(piston.gas().state())
            .map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || area > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid(
                "conducting film exceeds the finite body's exposed area",
            ));
        }
        let next = piston
            .thermal_candidate(
                input,
                body,
                controls,
                exchange,
                |gas, wall| {
                    let sample = self.resolve(gas, wall, || {
                        if (cancelled.borrow_mut())() {
                            Err(DepositionError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })?;
                    sample
                        .source
                        .admit_work(body.duration_s, exchange.source_tolerance_j)?;
                    Ok(sample)
                },
                || (cancelled.borrow_mut())(),
            )
            .map_err(|error| match error {
                ThermalPistonError::Interface(ConjugateGasError::Thermal(
                    DepositionError::Cancelled,
                )) => ThermalPistonError::Cancelled,
                other => other,
            })?;
        self.path
            .admit_state(next.thermal().state)
            .map_err(ThermalPistonError::Interface)?;
        if (cancelled.borrow_mut())() {
            return Err(ThermalPistonError::Cancelled);
        }
        Ok(next)
    }

    /// Use the original fixed-volume gas/body endpoint and wall-feedback owner.
    /// Its wall variable remains the BODY trace; each callback eliminates a
    /// separate physical surface through the original k(T)/film kernel.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = F::State>>(
        &self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<
                WallConvectionSample<ConjugateGasConvectionSample<F::Sample, Temperature>>,
                G,
            >,
        >,
        ConvectionError<WallConvectionError<ConjugateGasError>>,
    > {
        let source_error = |error| ConvectionError::Convection(WallConvectionError::Source(error));
        cx.checkpoint()
            .map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        self.path.admit_input(body).map_err(source_error)?;
        let area = self.area(environment.gas.state()).map_err(source_error)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || (area != 0.0 && area != body.surface_area_m2)
        {
            return Err(source_error(ConjugateGasError::Invalid(
                "isochoric conducting film requires full body coverage or explicit insulation",
            )));
        }
        if !environment.body.transport.radiation_is_disabled() {
            return Err(source_error(ConjugateGasError::Invalid(
                "conducting gas surfaces require disabled body-node radiation",
            )));
        }
        let next = environment
            .advance_with_wall_convection(cx, body, wall, |gas, body_temperature| {
                let sample = self.resolve_at_temperature(
                    gas,
                    body_temperature,
                    Temperature::new(body_temperature),
                    || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
                )?;
                sample
                    .source
                    .admit_work(body.duration_s, environment.exchange.source_tolerance_j)?;
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: if area == 0.0 {
                        0.0
                    } else {
                        sample.conductance_w_k / area
                    },
                    source: sample.source,
                })
            })
            .map_err(|error| match error {
                ConvectionError::Convection(WallConvectionError::Source(
                    ConjugateGasError::Thermal(DepositionError::Cancelled),
                )) => ConvectionError::Thermal(ThermalPistonError::Cancelled),
                other => other,
            })?;
        self.path.admit_state(next.state).map_err(source_error)?;
        cx.checkpoint()
            .map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        Ok(next)
    }
}

impl<'g, 'a, 'm, M: FilmGasSource> ConjugateGasSurface<'a, 'm, ConjugateFilmSelection<'g, M>> {
    pub(super) fn from_film(
        film: FilmGasSurface<'g, M>,
        flow: Option<MassFlowInput>,
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ReservoirLayerGeometry,
        config: LayerConductionConfig,
    ) -> Result<Self, ConjugateGasError> {
        Self::new(
            ConjugateFilmSelection { film, flow },
            layer,
            geometry,
            config,
        )
    }
}

impl<'g, M: FilmGasSource> FilmGasSurface<'g, M> {
    /// Bind the original pure-phase conducting path to this source film.
    /// The positive film area is checked against the physical face at each
    /// actual gas trial; native use also checks the exact body caloric object.
    pub fn with_conducting_layer<'a, 'm>(
        self,
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ReservoirLayerGeometry,
        config: LayerConductionConfig,
    ) -> Result<ConjugateFilmGasSurface<'g, 'a, 'm, M>, ConjugateGasError> {
        ConjugateFilmGasSurface::from_film(self, None, layer, geometry, config)
    }
}
