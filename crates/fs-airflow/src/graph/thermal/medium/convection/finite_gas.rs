//! Reciprocal convection between a finite material body and the actual gas cavity.
//!
//! The existing correlation owner supplies Nu and h from the CURRENT bulk gas
//! state supplied by the coupling owner: an energy midpoint for the piston or
//! an endpoint for the isochoric body. That owner solves one signed heat
//! transfer, body enthalpy (including latent heat), gas energy and any wall momentum.
//! No fixed-temperature reservoir, second gas law, or convection integrator.

/// Source-pinned film properties on the same finite gas/body heat residual.
pub mod film;

/// Multiple finite gases and convective patches on the existing circuit/thermal assembly.
pub mod network;
/// Buoyancy-driven vertical-plate heat exchange at actual gas/body states.
pub mod natural;
pub use natural::NaturalGasSurface;

use super::{ForcedConvectionSettings, MaterialConvectionError, evaluate_forced_from_material};
use fs_convection::{CorrelationId, HeatTransferCoefficient, NusseltEvaluation, ThermalDirection};
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_couple::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonCandidate, ThermalPistonError,
};
use fs_couple::finite_gas_thermal::{
    ConvectionCoefficient, ConvectionError, FiniteGasThermalEnvironment, FiniteGasThermalReport,
};
use fs_evidence::Evidence;
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::storage::IsochoricGasStorage;
use fs_qty::Length;

/// Explicit flow/geometry for a single exposed patch, not a stored coefficient.
/// Properties are evaluated at the gas BULK state supplied by the coupling
/// owner (piston energy midpoint or isochoric endpoint). Select a correlation
/// and temperature-contrast regime admitting that reference; no film-temperature
/// or stagnant-gas model is silently substituted. The prescribed speed is the
/// relative flow magnitude over this interval, not the translating wall velocity.
#[derive(Clone, Copy, Debug)]
pub struct ForcedGasSurface {
    pub area_m2: f64,
    pub characteristic_length: Length,
    pub speed_m_s: f64,
    pub correlation: CorrelationId,
    /// Direction is replaced by the actual hot/cold ordering when unequal.
    pub settings: ForcedConvectionSettings,
}

/// Original correlation evaluation and the SINGLE transport snapshot it used.
#[derive(Debug)]
pub struct GasConvectionEvaluation {
    pub transport: FluidTransportState,
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
}

/// Complete accepted interface sample; gas sources are not reduced to a tuple.
/// The body source identity remains in its canonical phase state/enthalpy chart.
#[derive(Debug)]
pub struct FiniteGasConvectionSample {
    pub gas: ResolvedThermochemicalGasState,
    pub body: EquilibriumPhaseState,
    pub area_m2: f64,
    /// None for an explicitly isolated zero-area patch or the natural law's
    /// exact thermal equilibrium. Neither case extrapolates a correlation.
    pub evaluation: Option<GasConvectionEvaluation>,
}

fn invalid(message: &str) -> MaterialConvectionError {
    MaterialError::Parameters { what: message.into() }.into()
}

impl ForcedGasSurface {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
            || !self.speed_m_s.is_finite() || self.speed_m_s < 0.0
        {
            return Err(invalid("finite gas convection requires nonnegative area/speed and positive length"));
        }
        Ok(())
    }

    /// Resolve the shared gas midpoint supplied by the ORIGINAL heat/motion
    /// residual. Pressure, composition, Cp, viscosity and conductivity all come
    /// from that same source-bearing state; there is no captured initial gas.
    /// This repeatable function is also usable by existing reactive, powered or
    /// multi-boundary callbacks. Geometry coverage belongs to their owner.
    ///
    /// G=h*A and Q=dt*G*(T_body-T_gas). This returns G and its sources only; the
    /// existing coupled solver owns Q and debits both sides exactly once.
    pub fn resolve(
        self, gas: &ResolvedThermochemicalGasState, body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FiniteGasConvectionSample>, MaterialConvectionError> {
        self.validate()?;
        let temperature = body.temperature_k();
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(invalid("finite gas convection needs a positive body temperature"));
        }
        if self.area_m2 == 0.0 {
            return Ok(ConductanceSample { conductance_w_k: 0.0,
                source: FiniteGasConvectionSample { gas: gas.clone(), body,
                    area_m2: 0.0, evaluation: None } });
        }
        let mut settings = self.settings;
        if temperature > gas.state().temperature {
            settings.direction = ThermalDirection::HeatingFluid;
        } else if temperature < gas.state().temperature {
            settings.direction = ThermalDirection::CoolingFluid;
        }
        let evaluated = evaluate_forced_from_material(self.correlation, gas,
            self.speed_m_s, self.characteristic_length, settings)?;
        let (transport, nusselt, coefficient) = evaluated.into_parts();
        let conductance = self.area_m2 * coefficient.value.value();
        if !conductance.is_finite() || conductance <= 0.0 {
            return Err(invalid("finite gas convective conductance overflowed or vanished"));
        }
        Ok(ConductanceSample { conductance_w_k: conductance,
            source: FiniteGasConvectionSample { gas: gas.clone(), body, area_m2: self.area_m2,
                evaluation: Some(GasConvectionEvaluation { transport, nusselt, coefficient }) } })
    }

    /// Couple a finite source-backed gas, inertial piston and finite enthalpy body.
    /// Returned endpoints must be published together; neither input is mutated.
    /// The patch must fit the declared body area and must not be counted in a
    /// second boundary law. External heat excludes this internal exchange.
    ///
    /// Uniform gas/body and negligible body pV/deformation work are the original
    /// owners' assumptions. Melting changes thermal state, NOT rigid wall inertia.
    /// Velocity is prescribed flow data: no jet momentum, fan energy, viscous
    /// heating, turbulence field, or boundary-layer resolution is inferred.
    /// Correlation-domain errors and insufficient temporal resolution refuse;
    /// no coefficient is clipped, extrapolated, or replaced with an air preset.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate(
        self, piston: &GasPiston, input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, exchange: HeatExchangeSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalPistonCandidate<FiniteGasConvectionSample>, ThermalPistonError<MaterialConvectionError>> {
        if cancelled() { return Err(ThermalPistonError::Cancelled); }
        self.validate().map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite() || body.surface_area_m2 <= 0.0
            || self.area_m2 > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid("convective patch exceeds the finite body's exposed area"));
        }
        piston.thermal_candidate(input, body, controls, exchange,
            |gas, wall| self.resolve(gas, wall), cancelled)
    }

    /// Reciprocal convection for an isochoric gas and the existing lumped body
    /// callback, including hereditary/deforming-solid consumers of that callback.
    /// The SAME candidate gas endpoint supplies storage, transport, h and the
    /// returned acoustic state. A data-only gas swap therefore changes the heat
    /// transfer law as well as thermodynamic storage; no piston is required.
    /// Ideal and Peng--Robinson parcels share this residual; nonideal density,
    /// heat capacity and pressure remain those of the real-gas endpoint.
    ///
    /// A positive patch must cover the whole exposed body boundary. Smearing a
    /// partial-patch h over that boundary would understate the local Biot number.
    /// Zero area explicitly disables convection. Radiation keeps its independent
    /// enclosure and the original body owner still admits conductivity and Biot.
    /// Bulk-endpoint property evaluation and prescribed relative speed retain the
    /// same correlation/domain assumptions as `resolve`; no film-temperature,
    /// flow-field or fan-power solution is inferred. The proposed heat selects
    /// Dittus-Boelter's heating/cooling exponent. Its zero-heat seed uses the
    /// smaller branch coefficient so an authored direction cannot cause a false
    /// Biot refusal before the physical transfer is known. Every nonzero trial
    /// uses its actual branch. Constant-heat-flux duct correlations remain
    /// incompatible with this isothermal-body model.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage>(
        self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
    ) -> Result<UniformEnthalpyStep<FiniteGasThermalReport<Option<GasConvectionEvaluation>, G>>,
        ConvectionError<MaterialConvectionError>>
    {
        cx.checkpoint().map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        self.validate().map_err(ConvectionError::Convection)?;
        if !body.surface_area_m2.is_finite() || body.surface_area_m2 <= 0.0
            || (self.area_m2 != 0.0 && self.area_m2 != body.surface_area_m2)
        {
            return Err(ConvectionError::Thermal(ThermalPistonError::Invalid(
                "isochoric convection must cover the complete body surface or be explicitly insulated",
            )));
        }
        if self.area_m2 > 0.0 && matches!(self.correlation,
            CorrelationId::CircularDuctLaminarChf | CorrelationId::RectangularDuctLaminarChf)
        {
            return Err(ConvectionError::Convection(invalid(
                "endpoint gas-only convection requires a wall-temperature-independent isothermal correlation",
            )));
        }
        environment.advance_with_convection_feedback(cx, body, |gas, heat_into_gas_j| {
            if self.area_m2 == 0.0 {
                return Ok(ConvectionCoefficient { coefficient_w_per_m2_k: 0.0, source: None });
            }
            let mut settings = self.settings;
            settings.direction = if heat_into_gas_j > 0.0 {
                ThermalDirection::HeatingFluid
            } else {
                ThermalDirection::CoolingFluid
            };
            let mut evaluated = evaluate_forced_from_material(self.correlation, gas,
                self.speed_m_s, self.characteristic_length, settings)?;
            if heat_into_gas_j == 0.0 && self.correlation == CorrelationId::DittusBoelter {
                settings.direction = ThermalDirection::HeatingFluid;
                let heating = evaluate_forced_from_material(self.correlation, gas,
                    self.speed_m_s, self.characteristic_length, settings)?;
                if heating.coefficient().value.value() < evaluated.coefficient().value.value() {
                    evaluated = heating;
                }
            }
            let (transport, nusselt, coefficient) = evaluated.into_parts();
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k: coefficient.value.value(),
                source: Some(GasConvectionEvaluation { transport, nusselt, coefficient }),
            })
        })
    }
}
