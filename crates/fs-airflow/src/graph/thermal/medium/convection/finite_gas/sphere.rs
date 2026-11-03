//! Source-pinned forced-sphere heat exchange on the existing finite heat solver.
//!
//! The same trial gas owns storage, pressure work and bulk transport. The
//! original source supplies fluid viscosity at the trial body temperature;
//! this second property query has no mass or energy of its own.

use super::film::FilmGasSource;
use super::*;
use crate::graph::thermal::medium::convection::evaluate_sphere_from_material;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// Exposed patch of an isolated isothermal sphere in a prescribed free stream.
/// Surface coverage belongs to the body owner. The sphere correlation includes
/// its conductive contribution once; do not add a parallel stagnant-gas link.
pub struct ForcedSphereGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    pub model: &'a M,
    pub area_m2: f64,
    pub diameter: Length,
    pub speed_m_s: f64,
}
impl<M: FilmGasSource> Copy for ForcedSphereGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for ForcedSphereGasSurface<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Actual coupled gas/body states and the separate fluid-at-surface query.
#[derive(Debug)]
pub struct SphereGasConvectionSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> {
    pub interface: FiniteGasConvectionSample<S, W>,
    /// None only for an explicitly insulated zero-area patch.
    pub surface_gas: Option<S>,
}

impl<M: FilmGasSource> ForcedSphereGasSurface<'_, M> {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        let diameter = self.diameter.value();
        if !self.area_m2.is_finite()
            || self.area_m2 < 0.0
            || !diameter.is_finite()
            || diameter <= 0.0
            || !self.speed_m_s.is_finite()
            || self.speed_m_s < 0.0
            || self.area_m2 > core::f64::consts::PI * diameter * diameter
        {
            return Err(invalid(
                "sphere convection needs positive diameter, nonnegative speed and an area no larger than pi D squared",
            ));
        }
        Ok(())
    }

    /// Resolve the original trial's bulk state and wall-temperature viscosity.
    /// Pressure, active composition and source claims stay pinned by the same
    /// ideal or explicit-branch Peng--Robinson material owner. `h` uses bulk k.
    ///
    /// The shared source range is enforced at every trial; ordinary gases with
    /// a hotter wall can have mu_bulk/mu_surface < 1 and are refused. This is a
    /// forced-flow sphere model, not a fallback for mixed convection, dense
    /// particle arrays, phase-changing boundary layers or a stagnant reservoir.
    /// Use this repeatable callback with existing powered/reactive/network
    /// heat owners, which retain the complete area and energy accounting.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<SphereGasConvectionSample<M::State>>, MaterialConvectionError>
    {
        self.resolve_at_temperature(gas, body.temperature_k(), body)
    }

    fn resolve_at_temperature<W>(
        self,
        gas: &M::State,
        wall: f64,
        body: W,
    ) -> Result<ConductanceSample<SphereGasConvectionSample<M::State, W>>, MaterialConvectionError>
    {
        self.validate()?;
        if !wall.is_finite() || wall <= 0.0 {
            return Err(invalid(
                "sphere convection needs a positive body temperature",
            ));
        }
        if self.area_m2 == 0.0 {
            return Ok(ConductanceSample {
                conductance_w_k: 0.0,
                source: SphereGasConvectionSample {
                    interface: FiniteGasConvectionSample {
                        gas: gas.clone(),
                        body,
                        area_m2: 0.0,
                        evaluation: None,
                    },
                    surface_gas: None,
                },
            });
        }
        let surface_gas = self
            .model
            .evaluate_property_temperature(gas, Temperature::new(wall))?;
        let evaluation =
            evaluate_sphere_from_material(gas, &surface_gas, self.speed_m_s, self.diameter)?;
        let (transport, nusselt, coefficient) = evaluation.into_parts();
        let conductance_w_k = self.area_m2 * coefficient.value.value();
        if !conductance_w_k.is_finite() || conductance_w_k <= 0.0 {
            return Err(invalid("sphere conductance overflowed or vanished"));
        }
        Ok(ConductanceSample {
            conductance_w_k,
            source: SphereGasConvectionSample {
                interface: FiniteGasConvectionSample {
                    gas: gas.clone(),
                    body,
                    area_m2: self.area_m2,
                    evaluation: Some(GasConvectionEvaluation {
                        transport,
                        nusselt,
                        coefficient,
                    }),
                },
                surface_gas: Some(surface_gas),
            },
        })
    }

    /// Advance gas energy, body enthalpy and piston motion together through the
    /// existing reciprocal heat residual. Both endpoints are returned or neither.
    /// The prescribed stream speed does not imply a fan, drag or kinetic-energy
    /// source. The complete nonlinear bracket must stay inside the source domain.
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
        ThermalPistonCandidate<SphereGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        self.validate().map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || self.area_m2 > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid(
                "sphere patch exceeds the finite body's exposed area",
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

    /// Resolve source-pinned surface viscosity at the coupled wall endpoint
    /// while the original fixed-volume gas retains its bulk transport and
    /// internal energy. The original wall-feedback/reciprocal-heat solve owns
    /// both thermal endpoints, radiation, body enthalpy and Biot admission.
    ///
    /// A positive patch must cover the complete exposed body surface. The
    /// entire heat and wall-temperature brackets must admit Whitaker's source
    /// bounds, including mu_bulk/mu_surface; no ratio or heat is clipped.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = M::State>>(
        self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<WallConvectionSample<SphereGasConvectionSample<M::State, Temperature>>, G>,
        >,
        ConvectionError<WallConvectionError<MaterialConvectionError>>,
    > {
        cx.checkpoint().map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        self.validate().map_err(|error| {
            ConvectionError::Convection(WallConvectionError::Source(error))
        })?;
        if !body.surface_area_m2.is_finite() || body.surface_area_m2 <= 0.0
            || (self.area_m2 != 0.0 && self.area_m2 != body.surface_area_m2)
        {
            return Err(ConvectionError::Thermal(ThermalPistonError::Invalid(
                "isochoric sphere convection must cover the complete body surface or be explicitly insulated",
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
