//! Finite original caloric contact followed by the native slab work/remap.

use super::*;
use fs_material::phase::exchange::pair::{
    FinitePhaseBody, FinitePhaseExchange, FinitePhaseExchangeControls, FinitePhaseExchangeError,
    exchange_with_checkpoint,
};

/// Simultaneous endpoint contact and native mechanical-work balance.
pub mod coupled;

/// Independent source-trajectory, represented-storage and native moving gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabReservoirHeatControls {
    pub heat: FreeSurfaceSlabHeatControls,
    pub exchange: FinitePhaseExchangeControls,
    /// Maximum |represented reservoir Delta H + heat into liquid| [J].
    pub maximum_reservoir_energy_defect_j: f64,
    /// Maximum |actual liquid Delta H + reservoir Delta H + native work| [J].
    pub maximum_total_energy_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceSlabReservoirHeatError {
    Native(FreeSurfaceSlabHeatError),
    Exchange(FinitePhaseExchangeError),
}
impl core::fmt::Display for FreeSurfaceSlabReservoirHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-reservoir thermal slab: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceSlabReservoirHeatError {}
impl From<FreeSurfaceSlabHeatError> for FreeSurfaceSlabReservoirHeatError {
    fn from(error: FreeSurfaceSlabHeatError) -> Self {
        Self::Native(error)
    }
}
impl From<Error> for FreeSurfaceSlabReservoirHeatError {
    fn from(error: Error) -> Self {
        Self::Native(error.into())
    }
}
impl From<FinitePhaseExchangeError> for FreeSurfaceSlabReservoirHeatError {
    fn from(error: FinitePhaseExchangeError) -> Self {
        Self::Exchange(error)
    }
}

/// One immutable proposal; neither participant advances before `advance` succeeds.
pub struct FreeSurfaceSlabReservoirHeatPlan<'a, 's, 'm, 'r> {
    native: FreeSurfaceSlabHeatPlan<'a, 's, 'm>,
    reservoir: FinitePhaseBody<'r>,
    exchange: FinitePhaseExchange,
    controls: FreeSurfaceSlabReservoirHeatControls,
    conductance_w_k: f64,
    duration_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabReservoirHeatReport {
    pub duration_s: f64,
    pub conductance_w_k: f64,
    pub heat_to_liquid_j: f64,
    pub reservoir_enthalpy_change_j: f64,
    pub reservoir_energy_defect_j: f64,
    pub total_energy_defect_j: f64,
    /// Final liquid T minus the provisional thermal-contact T [K]. This is a
    /// split-stage correction, not a bound on the time-discretization error.
    pub mechanical_temperature_correction_k: f64,
}

/// Both accepted endpoints. Reservoir mass/chart are unchanged; its `initial`
/// field is now the accepted state to supply to the next interval.
#[derive(Clone)]
pub struct FreeSurfaceSlabReservoirHeatStep<'s, 'm, 'r> {
    pub heat: FreeSurfaceSlabHeatStep<'s, 'm>,
    pub reservoir: FinitePhaseBody<'r>,
    /// Provisional contact-stage liquid state and entropy ONLY. Mechanical
    /// work subsequently changes the liquid; this is not coupled entropy.
    pub thermal_stage: FinitePhaseExchange,
    pub report: FreeSurfaceSlabReservoirHeatReport,
}

impl<'a, 's, 'm, 'r> FreeSurfaceSlabReservoirHeatPlan<'a, 's, 'm, 'r> {
    #[must_use]
    pub const fn native(&self) -> &FreeSurfaceSlabHeatPlan<'a, 's, 'm> {
        &self.native
    }

    /// Original finite-body heat flow. The liquid state is a source-split
    /// intermediate; only `advance` returns the mechanical-corrected endpoint.
    #[must_use]
    pub const fn thermal_stage(&self) -> &FinitePhaseExchange {
        &self.exchange
    }

    pub fn advance<'next>(
        &'next self,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabReservoirHeatStep<'next, 'next, 'r>, FreeSurfaceSlabReservoirHeatError>
    where
        'm: 'next,
    {
        let heat = self.native.advance(gate, cx)?;
        let actual = &heat.volume.report;
        let reservoir_change = Dd::from_f64(self.reservoir.mass_kg)
            * (Dd::from_f64(self.exchange.states[1].specific_enthalpy_j_kg())
                - Dd::from_f64(self.reservoir.initial.specific_enthalpy_j_kg()));
        let into_liquid = -self.exchange.heat_first_to_second_j;
        let reservoir_defect = value(reservoir_change + Dd::from_f64(into_liquid));
        let total_defect = value(
            Dd::from_f64(actual.requested_enthalpy_increment_j)
                + Dd::from_f64(actual.enthalpy_defect_j)
                + Dd::from_f64(actual.required_mechanical_work_j)
                + reservoir_change,
        );
        if !reservoir_defect.is_finite()
            || !total_defect.is_finite()
            || reservoir_defect.abs() > self.controls.maximum_reservoir_energy_defect_j
            || total_defect.abs() > self.controls.maximum_total_energy_defect_j
        {
            return Err(invalid(
                "represented reservoir/native liquid energy does not close",
            ));
        }
        let report = FreeSurfaceSlabReservoirHeatReport {
            duration_s: self.duration_s,
            conductance_w_k: self.conductance_w_k,
            heat_to_liquid_j: into_liquid,
            reservoir_enthalpy_change_j: value(reservoir_change),
            reservoir_energy_defect_j: reservoir_defect,
            total_energy_defect_j: total_defect,
            mechanical_temperature_correction_k: self.native.endpoint().state().temperature_k()
                - self.exchange.states[0].temperature_k(),
        };
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabReservoirHeatStep {
            heat,
            reservoir: FinitePhaseBody {
                initial: self.exchange.states[1],
                ..self.reservoir
            },
            thermal_stage: self.exchange,
            report,
        })
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Reciprocal finite-capacity heating/cooling of this original liquid mass.
    ///
    /// This is explicitly a FIRST-ORDER thermal-then-mechanical source split.
    /// The original finite pair integrates contact with fixed G over the native
    /// dt, debiting one body and crediting the other with ONE Q. The existing
    /// slab root then solves M*Delta h + W = Q, resolves the endpoint sources,
    /// and its native pressure/remap advances the interface and same material.
    /// It uses one native inertia preparation and one accepted pressure step.
    /// The temporary contact liquid state is not a published fluid endpoint.
    ///
    /// This does not integrate simultaneous heat/mechanical-power feedback.
    /// Timestep refinement is required to assess that splitting error. The
    /// pair's entropy belongs only to the thermal stage; no full coupled
    /// entropy/passivity claim follows from it. Projection and kinetic mapping
    /// losses retain their own native accounts and never become thermal heat.
    ///
    /// Both original fixed-pressure enthalpy charts, masses and source domains
    /// remain authoritative; ambient p*dV is already included and not added
    /// again. The finite reservoir is a caloric body, not a rigid sealed-vessel
    /// EOS. Its geometry or additional mechanical work requires its own owner.
    /// The caller must justify uniform temperature in both bodies and fixed
    /// total G (for example a fixed planar contact area). Smaller timesteps
    /// cannot establish that spatial-lumping assumption. The liquid remains
    /// pure; the reservoir may traverse its original equilibrium latent path.
    /// No external power, mass transfer, general ALE or viscous evolution is
    /// introduced. Refusal/cancellation publishes neither accepted endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_uniform_slab_reservoir_heat<'a, 'r>(
        &'a self,
        anchored_face: Face3,
        reservoir: FinitePhaseBody<'r>,
        conductance_w_k: f64,
        controls: FreeSurfaceSlabReservoirHeatControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabReservoirHeatPlan<'a, 's, 'm, 'r>, FreeSurfaceSlabReservoirHeatError>
    {
        checkpoint(gate, cx)?;
        if [
            conductance_w_k,
            controls.maximum_reservoir_energy_defect_j,
            controls.maximum_total_energy_defect_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "finite nonnegative contact conductance and energy gates required",
            ));
        }
        let prepared =
            self.prepare_uniform_slab_heat_state(anchored_face, controls.heat, gate, cx)?;
        let duration_s = prepared.dt;
        let exposure = value(Dd::from_f64(conductance_w_k) * Dd::from_f64(duration_s));
        if !exposure.is_finite() || (conductance_w_k > 0.0 && exposure == 0.0) {
            return Err(invalid("finite-reservoir exposure is unrepresentable"));
        }
        let exchange = exchange_with_checkpoint(
            FinitePhaseBody {
                curve: self.caloric,
                initial: prepared.initial,
                mass_kg: value(prepared.mass),
            },
            reservoir,
            exposure,
            controls.exchange,
            || gate.is_requested() || cx.is_cancel_requested(),
        )?;
        if exchange.states[0].phase() != SolidLiquidPhase::Liquid {
            return Err(invalid(
                "finite-reservoir thermal stage leaves the liquid phase",
            ));
        }
        if exchange.energy_residual_j[0].abs() > controls.heat.maximum_heat_defect_j
            || exchange.energy_residual_j[1].abs() > controls.maximum_reservoir_energy_defect_j
        {
            return Err(invalid(
                "original finite-body heat cannot be represented within its gates",
            ));
        }
        let native = prepared.with_supplied_heat(-exchange.heat_first_to_second_j, gate, cx)?;
        Ok(FreeSurfaceSlabReservoirHeatPlan {
            native,
            reservoir,
            exchange,
            controls,
            conductance_w_k,
            duration_s,
        })
    }
}

fn invalid(what: &'static str) -> FreeSurfaceSlabReservoirHeatError {
    FreeSurfaceSlabHeatError::Invalid(what).into()
}
