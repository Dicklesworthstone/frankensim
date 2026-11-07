//! Source-pinned turbulent-tube heat exchange on the existing finite heat solver.
//!
//! The same trial gas owns storage, pressure work and bulk transport. The
//! original source supplies fluid viscosity at the trial body temperature;
//! this second property query has no mass or energy of its own.

use super::film::FilmGasSource;
use super::*;
use crate::graph::thermal::medium::convection::evaluate_tube_from_material;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// Inside-wall patch of a smooth, straight circular tube with a prescribed
/// bulk stream speed and approximately uniform wall temperature.
///
/// The existing Sieder--Tate engineering closure uses bulk Re, Pr and k, with
/// only viscosity evaluated at the actual wall temperature. The caller must
/// establish fully developed, single-phase Newtonian flow. Its liquid-source
/// provenance does not establish empirical gas accuracy. Strong density, heat
/// capacity or conductivity variation, compressibility, buoyancy and phase
/// change in the fluid remain outside this frozen-bulk-property treatment.
/// Body phase change still belongs to the original enthalpy owner.
#[derive(Debug)]
pub struct ForcedTubeGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    pub model: &'a M,
    pub area_m2: f64,
    /// Actual inside circular diameter, not an arbitrary hydraulic length.
    pub diameter: Length,
    /// Physical heated length; splitting wall patches does not shorten it.
    pub heated_length: Length,
    pub speed_m_s: f64,
}
impl<M: FilmGasSource> Copy for ForcedTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for ForcedTubeGasSurface<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Actual coupled gas/body states and the separate fluid-at-surface query.
#[derive(Debug)]
pub struct TubeGasConvectionSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> {
    pub interface: FiniteGasConvectionSample<S, W>,
    /// None only for an explicitly insulated zero-area patch.
    pub surface_gas: Option<S>,
}

impl<M: FilmGasSource> ForcedTubeGasSurface<'_, M> {
    pub(super) fn validate(self) -> Result<(), MaterialConvectionError> {
        let diameter = self.diameter.value();
        let heated_length = self.heated_length.value();
        let length_ratio = heated_length / diameter;
        let inside_area = core::f64::consts::PI * diameter * heated_length;
        if !self.area_m2.is_finite()
            || self.area_m2 < 0.0
            || !diameter.is_finite()
            || diameter <= 0.0
            || !heated_length.is_finite()
            || heated_length <= 0.0
            || !length_ratio.is_finite()
            || length_ratio <= 0.0
            || !inside_area.is_finite()
            || inside_area <= 0.0
            || self.area_m2 > inside_area
            || !self.speed_m_s.is_finite()
            || self.speed_m_s < 0.0
        {
            return Err(invalid(
                "tube convection needs positive inside diameter and heated length, nonnegative speed and a patch no larger than pi D L",
            ));
        }
        Ok(())
    }

    /// Resolve the original trial's bulk state and wall-temperature viscosity.
    /// Pressure, active composition and source claims stay pinned by the same
    /// ideal or explicit-branch Peng--Robinson material owner. `h` uses bulk k.
    ///
    /// The shared card enforces its Re, Pr, heated L/D and viscosity-ratio
    /// admissions on every trial, in both heating and cooling directions.
    /// No ratio is clamped and no laminar or stagnant-flow fallback is supplied.
    /// Use this repeatable callback with existing powered/reactive/network
    /// heat owners, which retain the complete area and energy accounting.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<TubeGasConvectionSample<M::State>>, MaterialConvectionError> {
        self.resolve_at_temperature(gas, body.temperature_k(), body)
    }

    fn resolve_at_temperature<W>(
        self,
        gas: &M::State,
        wall: f64,
        body: W,
    ) -> Result<ConductanceSample<TubeGasConvectionSample<M::State, W>>, MaterialConvectionError>
    {
        self.validate()?;
        if !wall.is_finite() || wall <= 0.0 {
            return Err(invalid("tube convection needs a positive body temperature"));
        }
        if self.area_m2 == 0.0 {
            return Ok(ConductanceSample {
                conductance_w_k: 0.0,
                source: TubeGasConvectionSample {
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
        let evaluation = evaluate_tube_from_material(
            gas,
            &surface_gas,
            self.speed_m_s,
            self.diameter,
            self.heated_length.value() / self.diameter.value(),
        )?;
        let (transport, nusselt, coefficient) = evaluation.into_parts();
        let conductance_w_k = self.area_m2 * coefficient.value.value();
        if !conductance_w_k.is_finite() || conductance_w_k <= 0.0 {
            return Err(invalid("tube conductance overflowed or vanished"));
        }
        Ok(ConductanceSample {
            conductance_w_k,
            source: TubeGasConvectionSample {
                interface: FiniteGasConvectionSample {
                    gas: gas.clone(),
                    body,
                    area_m2: self.area_m2,
                    evaluation: Some(GasConvectionEvaluation {
                        transport,
                        nusselt,
                        coefficient,
                        section: None,
                        tube_section: None,
                        rectangular_section: None,
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
        ThermalPistonCandidate<TubeGasConvectionSample<M::State>, G>,
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
                "tube patch exceeds the finite body's exposed area",
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
    /// entire heat and wall-temperature brackets must admit Sieder--Tate's source
    /// bounds, including mu_bulk/mu_surface; no ratio or heat is clipped.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = M::State>>(
        self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<
                WallConvectionSample<TubeGasConvectionSample<M::State, Temperature>>,
                G,
            >,
        >,
        ConvectionError<WallConvectionError<MaterialConvectionError>>,
    > {
        cx.checkpoint()
            .map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        self.validate()
            .map_err(|error| ConvectionError::Convection(WallConvectionError::Source(error)))?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || (self.area_m2 != 0.0 && self.area_m2 != body.surface_area_m2)
        {
            return Err(ConvectionError::Thermal(ThermalPistonError::Invalid(
                "isochoric tube convection must cover the complete body surface or be explicitly insulated",
            )));
        }
        environment.advance_with_wall_convection(cx, body, wall, |gas, temperature_k| {
            let sample =
                self.resolve_at_temperature(gas, temperature_k, Temperature::new(temperature_k))?;
            let coefficient_w_per_m2_k = sample
                .source
                .interface
                .evaluation
                .as_ref()
                .map_or(0.0, |value| value.coefficient.value.value());
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k,
                source: sample.source,
            })
        })
    }
}
