//! Source-pinned circular-tube heat exchange on the existing finite heat solver.
//!
//! The same trial gas owns storage, pressure work and bulk transport. The
//! original source supplies fluid viscosity at the trial body temperature;
//! this second property query has no mass or energy of its own.

use super::conjugate::{ConjugateGasFilm, sealed as conducting_sealed};
use super::film::FilmGasSource;
use super::mass_flow::MassFlowInput;
use super::*;
use crate::graph::thermal::medium::convection::{
    evaluate_laminar_tube_from_material, evaluate_tube_from_material,
};
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

/// Explicit whole-tube laminar Sieder--Tate mean on the SAME finite heat owners.
///
/// Construct with [`ForcedTubeGasSurface::laminar`]. The entire positive
/// heated area must be `pi D L`, with an approximately uniform wall temperature.
/// The original homogeneous trial gas represents the caller's arithmetic-mean
/// inlet/outlet bulk temperature required by the source relation. This adapter
/// does not resolve axial advection, inlet/outlet temperatures or a local profile;
/// splitting the tube cannot restart or shorten its physical heated length.
///
/// Original laminar, Prandtl, wall-viscosity and thermal-entry admissions remain
/// active on every trial. The empirical liquid-source relation does not establish
/// measured gas accuracy; no fully developed or turbulent fallback is introduced.
#[derive(Debug)]
pub struct LaminarTubeGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    surface: ForcedTubeGasSurface<'a, M>,
}
impl<M: FilmGasSource> Copy for LaminarTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for LaminarTubeGasSurface<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Original Sieder--Tate tube law with an explicit stream/recirculation rate.
///
/// Construct from the turbulent surface or its explicit `.laminar()` selection.
/// The actual circular flow area must equal `pi D^2/4`; the heated length and
/// laminar whole-wall admission remain unchanged. Each nonlinear gas trial
/// supplies original EOS density for velocity, while its separate same-source
/// wall query still supplies only viscosity. This adds no open mass/enthalpy
/// flux, fan work, axial solution or thermal integrator.
#[derive(Debug)]
pub struct PrescribedMassFlowTubeGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    surface: ForcedTubeGasSurface<'a, M>,
    law: TubeLaw,
    flow: MassFlowInput,
}
impl<M: FilmGasSource> Copy for PrescribedMassFlowTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for PrescribedMassFlowTubeGasSurface<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

#[derive(Clone, Copy, Debug)]
enum TubeLaw {
    Turbulent,
    Laminar,
}

/// Actual coupled gas/body states and the separate fluid-at-surface query.
#[derive(Debug)]
pub struct TubeGasConvectionSample<S = ResolvedThermochemicalGasState, W = EquilibriumPhaseState> {
    pub interface: FiniteGasConvectionSample<S, W>,
    /// None only for an explicitly insulated zero-area patch.
    pub surface_gas: Option<S>,
    /// Velocity from the ORIGINAL bulk trial when mass flow was selected.
    pub prescribed_mass_flow: Option<GasMassFlowState>,
}

impl<'a, M: FilmGasSource> ForcedTubeGasSurface<'a, M> {
    /// Explicitly select the original laminar mean over the complete heated
    /// tube. The native finite storage/heat solve and original source stay fixed.
    pub const fn laminar(self) -> LaminarTubeGasSurface<'a, M> {
        LaminarTubeGasSurface { surface: self }
    }

    /// Replace authored speed with positive mass flow through the same tube.
    pub fn with_prescribed_mass_flow(
        self,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<PrescribedMassFlowTubeGasSurface<'a, M>, MaterialConvectionError> {
        self.select_mass_flow(TubeLaw::Turbulent, mass_flow_kg_s, flow_area_m2)
    }

    fn select_mass_flow(
        self,
        law: TubeLaw,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<PrescribedMassFlowTubeGasSurface<'a, M>, MaterialConvectionError> {
        let flow = MassFlowInput::new(mass_flow_kg_s, flow_area_m2)?;
        flow.admit_circular_area(self.diameter.value())?;
        Ok(PrescribedMassFlowTubeGasSurface {
            surface: self,
            law,
            flow,
        })
    }

    fn bind_flow(
        &mut self,
        flow: Option<MassFlowInput>,
        gas: &M::State,
    ) -> Result<Option<GasMassFlowState>, MaterialConvectionError> {
        match flow {
            Some(flow) => {
                let state = flow.resolve(gas)?;
                self.speed_m_s = state.speed_m_s;
                Ok(Some(state))
            }
            None => Ok(None),
        }
    }

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
            || self.area_m2 - inside_area > 64.0 * f64::EPSILON * inside_area
            || !self.speed_m_s.is_finite()
            || self.speed_m_s < 0.0
        {
            return Err(invalid(
                "tube convection needs positive inside diameter and heated length, nonnegative speed and a patch no larger than pi D L",
            ));
        }
        Ok(())
    }

    fn validate_for(self, law: TubeLaw) -> Result<(), MaterialConvectionError> {
        self.validate()?;
        let inside_area =
            core::f64::consts::PI * self.diameter.value() * self.heated_length.value();
        if matches!(law, TubeLaw::Laminar)
            && self.area_m2 != 0.0
            && (self.area_m2 - inside_area).abs() > 64.0 * f64::EPSILON * inside_area
        {
            return Err(invalid(
                "laminar Sieder-Tate mean requires the complete physical heated area pi D L",
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
        self.resolve_at_temperature(TubeLaw::Turbulent, None, gas, body.temperature_k(), body)
    }

    fn resolve_at_temperature<W>(
        mut self,
        law: TubeLaw,
        flow: Option<MassFlowInput>,
        gas: &M::State,
        wall: f64,
        body: W,
    ) -> Result<ConductanceSample<TubeGasConvectionSample<M::State, W>>, MaterialConvectionError>
    {
        let prescribed_mass_flow = self.bind_flow(flow, gas)?;
        self.validate_for(law)?;
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
                    prescribed_mass_flow,
                },
            });
        }
        let surface_gas = self
            .model
            .evaluate_property_temperature(gas, Temperature::new(wall))?;
        let length_ratio = self.heated_length.value() / self.diameter.value();
        let evaluation = match law {
            TubeLaw::Turbulent => evaluate_tube_from_material(
                gas,
                &surface_gas,
                self.speed_m_s,
                self.diameter,
                length_ratio,
            )?,
            TubeLaw::Laminar => evaluate_laminar_tube_from_material(
                gas,
                &surface_gas,
                self.speed_m_s,
                self.diameter,
                length_ratio,
            )?,
        };
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
                prescribed_mass_flow,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<TubeGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        self.finite_body_candidate_with(
            TubeLaw::Turbulent,
            None,
            piston,
            input,
            body,
            controls,
            exchange,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn finite_body_candidate_with<G: PistonGas<State = M::State>>(
        mut self,
        law: TubeLaw,
        flow: Option<MassFlowInput>,
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
        self.bind_flow(flow, piston.gas().state())
            .map_err(ThermalPistonError::Interface)?;
        self.validate_for(law)
            .map_err(ThermalPistonError::Interface)?;
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
            |gas, wall| self.resolve_at_temperature(law, flow, gas, wall.temperature_k(), wall),
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
        self.isochoric_body_candidate_with(TubeLaw::Turbulent, None, cx, environment, body, wall)
    }

    fn isochoric_body_candidate_with<G: IsochoricGasStorage<State = M::State>>(
        mut self,
        law: TubeLaw,
        flow: Option<MassFlowInput>,
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
        self.bind_flow(flow, environment.gas.state())
            .map_err(|error| ConvectionError::Convection(WallConvectionError::Source(error)))?;
        self.validate_for(law)
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
            let sample = self.resolve_at_temperature(
                law,
                flow,
                gas,
                temperature_k,
                Temperature::new(temperature_k),
            )?;
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

impl<'a, M: FilmGasSource> LaminarTubeGasSurface<'a, M> {
    /// Keep the explicit laminar whole-tube mean while prescribing stream mass
    /// flow through its actual circular cross section.
    pub fn with_prescribed_mass_flow(
        self,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<PrescribedMassFlowTubeGasSurface<'a, M>, MaterialConvectionError> {
        self.surface
            .select_mass_flow(TubeLaw::Laminar, mass_flow_kg_s, flow_area_m2)
    }

    /// Evaluate the original whole-length mean using this trial's bulk state
    /// and same-fluid wall viscosity. No finite-section evidence is fabricated.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<TubeGasConvectionSample<M::State>>, MaterialConvectionError> {
        self.surface
            .resolve_at_temperature(TubeLaw::Laminar, None, gas, body.temperature_k(), body)
    }

    /// Advance the same finite gas, body enthalpy and piston momentum residual.
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
        ThermalPistonCandidate<TubeGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        self.surface.finite_body_candidate_with(
            TubeLaw::Laminar,
            None,
            piston,
            input,
            body,
            controls,
            exchange,
            cancelled,
        )
    }

    /// Advance the original fixed-volume gas and body wall-feedback solve.
    /// A positive whole-tube area must also cover the complete exposed body.
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
        self.surface.isochoric_body_candidate_with(
            TubeLaw::Laminar,
            None,
            cx,
            environment,
            body,
            wall,
        )
    }
}

impl<M: FilmGasSource> PrescribedMassFlowTubeGasSurface<'_, M> {
    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.flow.mass_flow_kg_s()
    }
    pub const fn flow_area_m2(&self) -> f64 {
        self.flow.flow_area_m2()
    }

    /// Resolve actual bulk velocity and original same-fluid wall viscosity.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<TubeGasConvectionSample<M::State>>, MaterialConvectionError> {
        self.surface.resolve_at_temperature(
            self.law,
            Some(self.flow),
            gas,
            body.temperature_k(),
            body,
        )
    }

    /// Use the same finite gas/body/piston residual with per-trial velocity.
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
        ThermalPistonCandidate<TubeGasConvectionSample<M::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        self.surface.finite_body_candidate_with(
            self.law,
            Some(self.flow),
            piston,
            input,
            body,
            controls,
            exchange,
            cancelled,
        )
    }

    /// Use the same fixed-volume gas/body and source-wall feedback residual.
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
        self.surface.isochoric_body_candidate_with(
            self.law,
            Some(self.flow),
            cx,
            environment,
            body,
            wall,
        )
    }
}

impl<M: FilmGasSource> conducting_sealed::Film for ForcedTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> ConjugateGasFilm for ForcedTubeGasSurface<'_, M> {
    type State = M::State;
    type Sample = TubeGasConvectionSample<M::State, Temperature>;

    fn admitted_area_m2(&self, _gas: &Self::State) -> Result<f64, MaterialConvectionError> {
        self.validate_for(TubeLaw::Turbulent)?;
        Ok(self.area_m2)
    }

    fn resolve_film_at_surface(
        &self,
        gas: &Self::State,
        surface: Temperature,
    ) -> Result<ConductanceSample<Self::Sample>, MaterialConvectionError> {
        self.resolve_at_temperature(TubeLaw::Turbulent, None, gas, surface.value(), surface)
    }
}

impl<M: FilmGasSource> conducting_sealed::Film for LaminarTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> ConjugateGasFilm for LaminarTubeGasSurface<'_, M> {
    type State = M::State;
    type Sample = TubeGasConvectionSample<M::State, Temperature>;

    fn admitted_area_m2(&self, _gas: &Self::State) -> Result<f64, MaterialConvectionError> {
        self.surface.validate_for(TubeLaw::Laminar)?;
        Ok(self.surface.area_m2)
    }

    fn resolve_film_at_surface(
        &self,
        gas: &Self::State,
        surface: Temperature,
    ) -> Result<ConductanceSample<Self::Sample>, MaterialConvectionError> {
        self.surface
            .resolve_at_temperature(TubeLaw::Laminar, None, gas, surface.value(), surface)
    }
}

impl<M: FilmGasSource> conducting_sealed::Film for PrescribedMassFlowTubeGasSurface<'_, M> {}
impl<M: FilmGasSource> ConjugateGasFilm for PrescribedMassFlowTubeGasSurface<'_, M> {
    type State = M::State;
    type Sample = TubeGasConvectionSample<M::State, Temperature>;

    fn admitted_area_m2(&self, gas: &Self::State) -> Result<f64, MaterialConvectionError> {
        let mut surface = self.surface;
        surface.bind_flow(Some(self.flow), gas)?;
        surface.validate_for(self.law)?;
        Ok(surface.area_m2)
    }

    fn resolve_film_at_surface(
        &self,
        gas: &Self::State,
        surface: Temperature,
    ) -> Result<ConductanceSample<Self::Sample>, MaterialConvectionError> {
        self.surface.resolve_at_temperature(
            self.law,
            Some(self.flow),
            gas,
            surface.value(),
            surface,
        )
    }
}
