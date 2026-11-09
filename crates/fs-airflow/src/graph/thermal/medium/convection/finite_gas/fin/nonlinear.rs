//! Original material and gas sources on a spatial nonlinear fin balance.
//!
//! Every axial edge uses the original exact mean of k(T), including crossed
//! source knots. Each surface node evaluates the selected gas film at its own
//! temperature. The film's authored length/flow regime remains a local closure
//! assumption: this does not solve the external boundary layer or fin-array
//! interference. The base body alone owns thermal storage; fins are steady.

use super::*;
use fs_convection::fin::nonlinear::{
    self as kernel, FinProperties, NonlinearFinControls, NonlinearFinError, NonlinearFinResponse,
    NonlinearFinTip,
};
use std::cell::RefCell;

/// Original source refusal or the spatial fin solver's bounded failure.
#[derive(Debug)]
pub enum NonlinearSourceFinError {
    Source(SourceFinError),
    Solve(NonlinearFinError<SourceFinError>),
}
impl core::fmt::Display for NonlinearSourceFinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "distributed source fin convection: {self:?}")
    }
}
impl std::error::Error for NonlinearSourceFinError {}
impl From<SourceFinError> for NonlinearSourceFinError {
    fn from(error: SourceFinError) -> Self {
        Self::Source(error)
    }
}

/// Actual accepted axial profile and every original local gas-film payload.
/// The node-zero side sample also supplies the bare base film, without another
/// source query. Samples retain their original authored area as diagnostics;
/// the numerical fin owner applies the local perimeter/tip area exactly once.
#[derive(Debug)]
pub struct NonlinearSourceFinSample<S> {
    pub conductivity_source: Arc<PropertyLineSlice>,
    pub geometry: StraightFinArray,
    pub base_temperature_k: f64,
    pub bulk_temperature_k: f64,
    pub fin: NonlinearFinResponse<S>,
    pub bare_base_conductance_w_k: f64,
    pub total_conductance_w_k: f64,
    /// Bare base convection plus the heat drawn at all fin roots [W].
    pub applied_heat_into_gas_w: f64,
    /// Root heat minus integrated fin side/tip heat, summed over the array [W].
    pub array_heat_residual_w: f64,
    /// Sum of absolute local spatial residuals, summed over the array [W].
    pub array_absolute_heat_defect_w: f64,
}

/// Source-pinned steady fins with k(T) and h(T_surface,T_bulk) resolved locally.
///
/// The original all-solid fin material may differ from the base body's material.
/// `controls.absolute_heat_tolerance_w` applies to the COMPLETE repeated array;
/// it is divided by the fin count before the single representative fin solve.
/// The caller owns spatial refinement and the local film model's applicability.
/// Transverse Biot admission is checked at the evaluated local surface states,
/// not certified over unsampled temperatures or an unresolved boundary layer.
pub struct NonlinearSourceFinGasSurface<'a, 'm, F: ConjugateGasFilm> {
    film: F,
    material: &'a ResolvedConductionLayer<'m>,
    geometry: StraightFinArray,
    controls: NonlinearFinControls,
}

impl<'a, 'm, F: ConjugateGasFilm> NonlinearSourceFinGasSurface<'a, 'm, F> {
    pub fn new(
        film: F,
        material: &'a ResolvedConductionLayer<'m>,
        geometry: StraightFinArray,
        controls: NonlinearFinControls,
    ) -> Result<Self, NonlinearSourceFinError> {
        controls.validate().map_err(SourceFinError::Invalid)?;
        if material.phase() != SolidLiquidPhase::Solid {
            return Err(SourceFinError::Invalid(
                "a fixed fin requires an original all-solid conductivity branch",
            )
            .into());
        }
        Ok(Self {
            film,
            material,
            geometry,
            controls,
        })
    }

    pub const fn film(&self) -> &F {
        &self.film
    }
    pub const fn material(&self) -> &ResolvedConductionLayer<'m> {
        self.material
    }
    pub const fn geometry(&self) -> StraightFinArray {
        self.geometry
    }
    pub const fn controls(&self) -> NonlinearFinControls {
        self.controls
    }

    /// Resolve one immutable bulk/base trial. Rejected iterations consume the
    /// original solver's source/work allowance but publish no physical history.
    pub fn resolve(
        &self,
        gas: &F::State,
        base: Temperature,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<NonlinearSourceFinSample<F::Sample>>, NonlinearSourceFinError>
    {
        self.resolve_with_controls(gas, base, self.controls, checkpoint)
    }

    fn resolve_with_controls(
        &self,
        gas: &F::State,
        base: Temperature,
        mut controls: NonlinearFinControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<NonlinearSourceFinSample<F::Sample>>, NonlinearSourceFinError>
    {
        checkpoint().map_err(SourceFinError::Thermal)?;
        let geometry = self.geometry;
        let area = self
            .film
            .admitted_area_m2(gas)
            .map_err(SourceFinError::Film)?;
        if area != geometry.wetted_area_m2() {
            return Err(SourceFinError::Invalid(
                "fin film must cover the actual bare-base plus wetted-fin area exactly",
            )
            .into());
        }
        let count = f64::from(geometry.count());
        controls.absolute_heat_tolerance_w /= count;
        controls.validate().map_err(SourceFinError::Invalid)?;
        let curve = self.material.conductivity_curve();
        let stopped = RefCell::new(None);
        let response = kernel::solve(
            geometry.fin(),
            base.value(),
            gas.temperature_k(),
            match geometry.tip() {
                FinTipExposure::Adiabatic => NonlinearFinTip::Adiabatic,
                FinTipExposure::Exposed => NonlinearFinTip::Convection,
            },
            controls,
            |point| {
                let conductivity_w_m_k = curve.conductivity_w_m_k(point.temperature_k)?;
                let film = self
                    .film
                    .resolve_film_at_surface(gas, Temperature::new(point.temperature_k))?;
                Ok(FinProperties {
                    conductivity_w_m_k,
                    coefficient_w_m2_k: film.conductance_w_k / area,
                    source: film.source,
                })
            },
            |a, b| {
                curve
                    .mean_conductivity_w_m_k(a, b)
                    .map_err(SourceFinError::Thermal)
            },
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    *stopped.borrow_mut() = Some(error);
                    true
                }
            },
        );
        // Preserve the original checkpoint's typed refusal, including a caller
        // error other than cancellation; the numerical core only sees a stop.
        if let Some(error) = stopped.into_inner() {
            return Err(SourceFinError::Thermal(error).into());
        }
        let fin = response.map_err(NonlinearSourceFinError::Solve)?;
        let first = fin
            .nodes
            .first()
            .ok_or(SourceFinError::Fin(FinError::Arithmetic))?;
        let bare = first.properties.coefficient_w_m2_k * geometry.bare_base_area_m2();
        let conductance = bare + count * fin.effective_conductance_w_k;
        let bulk = gas.temperature_k();
        let heat = conductance * (base.value() - bulk);
        let direct_heat = bare * (base.value() - bulk) + count * fin.base_heat_w;
        let residual = count * fin.energy_residual_w;
        let defect = count * fin.sum_absolute_node_residual_w;
        let rounding = 32.0 * f64::EPSILON * heat.abs().max(direct_heat.abs());
        if ![bare, conductance, heat, direct_heat, residual, defect]
            .into_iter()
            .all(f64::is_finite)
            || bare < 0.0
            || conductance < 0.0
            || defect < 0.0
            || (heat - direct_heat).abs() > rounding
            || (conductance > 0.0 && base.value() != bulk && heat == 0.0)
        {
            return Err(SourceFinError::Fin(FinError::Arithmetic).into());
        }
        checkpoint().map_err(SourceFinError::Thermal)?;
        Ok(ConductanceSample {
            conductance_w_k: conductance,
            source: NonlinearSourceFinSample {
                conductivity_source: Arc::clone(self.material.source_line()),
                geometry,
                base_temperature_k: base.value(),
                bulk_temperature_k: bulk,
                fin,
                bare_base_conductance_w_k: bare,
                total_conductance_w_k: conductance,
                applied_heat_into_gas_w: heat,
                array_heat_residual_w: residual,
                array_absolute_heat_defect_w: defect,
            },
        })
    }

    /// One effective coefficient on the original native BASE area. Local h*A
    /// samples are evidence and must not be independently deposited as heat.
    pub fn resolve_isochoric_coefficient(
        &self,
        cx: &Cx<'_>,
        gas: &F::State,
        body_surface_area_m2: f64,
        base: Temperature,
    ) -> Result<ConvectionCoefficient<NonlinearSourceFinSample<F::Sample>>, NonlinearSourceFinError>
    {
        self.coefficient_with_controls(cx, gas, body_surface_area_m2, base, self.controls)
    }

    fn coefficient_with_controls(
        &self,
        cx: &Cx<'_>,
        gas: &F::State,
        area: f64,
        base: Temperature,
        controls: NonlinearFinControls,
    ) -> Result<ConvectionCoefficient<NonlinearSourceFinSample<F::Sample>>, NonlinearSourceFinError>
    {
        if area != self.geometry.base_area_m2() {
            return Err(SourceFinError::Invalid(
                "native fin body must cover the original base surface exactly",
            )
            .into());
        }
        let sample = self.resolve_with_controls(gas, base, controls, || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        })?;
        let coefficient = sample.conductance_w_k / area;
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err(SourceFinError::Fin(FinError::Arithmetic).into());
        }
        Ok(ConvectionCoefficient {
            coefficient_w_per_m2_k: coefficient,
            source: sample.source,
        })
    }

    /// Use the original reciprocal finite-gas/body solve, evaluating the entire
    /// nonlinear fin at every actual wall/gas trial. Array spatial residuals
    /// receive at most the interval's existing source-energy allowance divided
    /// by its actual duration. A relative tolerance cannot widen that allowance.
    ///
    /// The native body owns only its declared inventory. No fin mass, storage,
    /// radiation, contact resistance, geometry evolution or fan work is inferred.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = F::State>>(
        &self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<WallConvectionSample<NonlinearSourceFinSample<F::Sample>>, G>,
        >,
        ConvectionError<WallConvectionError<NonlinearSourceFinError>>,
    > {
        let source_error = |e| match e {
            NonlinearSourceFinError::Source(SourceFinError::Thermal(
                DepositionError::Cancelled,
            )) => ConvectionError::Thermal(ThermalPistonError::Cancelled),
            other => ConvectionError::Convection(WallConvectionError::Source(other)),
        };
        cx.checkpoint()
            .map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        if !environment.body.transport.radiation_is_disabled() {
            return Err(source_error(
                SourceFinError::Invalid("quasi-steady fins require disabled body-node radiation")
                    .into(),
            ));
        }
        let allowance = environment.exchange.source_tolerance_j / body.duration_s;
        if !body.duration_s.is_finite()
            || body.duration_s <= 0.0
            || !allowance.is_finite()
            || allowance <= 0.0
        {
            return Err(source_error(
                SourceFinError::Invalid(
                    "nonlinear fin feedback needs positive duration and source-energy allowance",
                )
                .into(),
            ));
        }
        let mut controls = self.controls;
        controls.absolute_heat_tolerance_w = controls.absolute_heat_tolerance_w.min(allowance);
        controls.relative_heat_tolerance = 0.0;
        self.coefficient_with_controls(
            cx,
            environment.gas.state(),
            body.surface_area_m2,
            Temperature::new(body.initial.temperature_k()),
            controls,
        )
        .map_err(source_error)?;
        environment
            .advance_with_wall_convection(cx, body, wall, |gas, t| {
                self.coefficient_with_controls(
                    cx,
                    gas,
                    body.surface_area_m2,
                    Temperature::new(t),
                    controls,
                )
            })
            .map_err(|error| match error {
                ConvectionError::Convection(WallConvectionError::Source(
                    NonlinearSourceFinError::Source(SourceFinError::Thermal(
                        DepositionError::Cancelled,
                    )),
                )) => ConvectionError::Thermal(ThermalPistonError::Cancelled),
                other => other,
            })
    }
}

impl<'g, M: FilmGasSource> FilmGasSurface<'g, M> {
    /// Resolve the original film locally along source-dependent steady fins.
    pub fn with_nonlinear_straight_fins<'a, 'm>(
        self,
        material: &'a ResolvedConductionLayer<'m>,
        geometry: StraightFinArray,
        controls: NonlinearFinControls,
    ) -> Result<
        NonlinearSourceFinGasSurface<'a, 'm, ConjugateFilmSelection<'g, M>>,
        NonlinearSourceFinError,
    > {
        NonlinearSourceFinGasSurface::new(
            ConjugateFilmSelection::new(self, None),
            material,
            geometry,
            controls,
        )
    }
}
