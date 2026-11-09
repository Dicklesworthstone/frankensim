//! Original finite caloric exchange coupled to native slab expansion work.

use super::*;
use fs_material::phase::exchange::{PhaseDrivenHeatExchange, PhaseHeatExchangeError};

/// Authored uniform heating and one fixed-temperature reservoir. Conductance
/// is the complete physical G [W/K], held fixed over this native interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabBath {
    pub reservoir_temperature_k: f64,
    pub conductance_w_k: f64,
    pub applied_power_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceSlabBathHeatError {
    Native(FreeSurfaceSlabHeatError),
    Exchange(PhaseHeatExchangeError),
}
impl core::fmt::Display for FreeSurfaceSlabBathHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-bath thermal slab: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceSlabBathHeatError {}
impl From<FreeSurfaceSlabHeatError> for FreeSurfaceSlabBathHeatError {
    fn from(error: FreeSurfaceSlabHeatError) -> Self {
        Self::Native(error)
    }
}
impl From<Error> for FreeSurfaceSlabBathHeatError {
    fn from(error: Error) -> Self {
        Self::Native(error.into())
    }
}
impl From<PhaseHeatExchangeError> for FreeSurfaceSlabBathHeatError {
    fn from(error: PhaseHeatExchangeError) -> Self {
        Self::Exchange(error)
    }
}

/// One admitted source/pressure-work root, borrowing the unchanged native slab.
pub struct FreeSurfaceSlabBathHeatPlan<'a, 's, 'm> {
    native: FreeSurfaceSlabHeatPlan<'a, 's, 'm>,
    source: FreeSurfaceSlabBath,
    exchange: PhaseDrivenHeatExchange,
    duration_s: f64,
    applied_heat: Dd,
    reservoir_heat: Dd,
}

/// Accepted moving native liquid and independently measured external heat balance.
#[derive(Clone)]
pub struct FreeSurfaceSlabBathHeatStep<'s, 'm> {
    pub heat: FreeSurfaceSlabHeatStep<'s, 'm>,
    pub source: FreeSurfaceSlabBath,
    pub duration_s: f64,
    pub mean_liquid_temperature_k: f64,
    pub applied_heat_j: f64,
    /// G*dt*(T_bath - mean T), signed into the same conserved liquid mass.
    pub reservoir_heat_j: f64,
    /// Actual native Delta H + correction work - applied heat - reservoir heat.
    pub heat_balance_defect_j: f64,
}

impl<'a, 's, 'm> FreeSurfaceSlabBathHeatPlan<'a, 's, 'm> {
    #[must_use]
    pub const fn native(&self) -> &FreeSurfaceSlabHeatPlan<'a, 's, 'm> {
        &self.native
    }

    /// Original caloric flow at the accepted root trial. Its continuous heat
    /// and native-state rounding remain distinct; the moving owner audits the
    /// actual represented endpoint against the physical joule allowance.
    #[must_use]
    pub const fn thermal_exchange(&self) -> &PhaseDrivenHeatExchange {
        &self.exchange
    }

    /// Execute the original native pressure/remap once. No source debit,
    /// liquid clock or thermal state is published on refusal or cancellation.
    pub fn advance<'next>(
        &'next self,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabBathHeatStep<'next, 'next>, FreeSurfaceSlabBathHeatError>
    where
        'm: 'next,
    {
        let heat = self.native.advance(gate, cx)?;
        let actual = &heat.volume.report;
        let defect = value(
            Dd::from_f64(actual.requested_enthalpy_increment_j)
                + Dd::from_f64(actual.enthalpy_defect_j)
                + Dd::from_f64(actual.required_mechanical_work_j)
                - self.applied_heat
                - self.reservoir_heat,
        );
        if !defect.is_finite() || defect.abs() > self.native.controls.maximum_heat_defect_j {
            return Err(invalid(
                "actual moving slab exceeds its integrated bath/source heat budget",
            ));
        }
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabBathHeatStep {
            heat,
            source: self.source,
            duration_s: self.duration_s,
            mean_liquid_temperature_k: self.exchange.mean_temperature_k,
            applied_heat_j: value(self.applied_heat),
            reservoir_heat_j: value(self.reservoir_heat),
            heat_balance_defect_j: defect,
        })
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Heat/cool the same finite liquid mass from a bath and an applied power,
    /// paying the ORIGINAL native pressure-correction work during the interval.
    ///
    /// The explicitly declared source split holds that mechanical power at
    /// W/dt while the original caloric owner solves
    /// m*dh/dt=P+G*(T_bath-T(h))-W/dt. For each endpoint trial, W=B*u*(u-u0)
    /// comes from the same measured native face inertia and source-volume
    /// velocity as `prepare_uniform_slab_heat`. Its bounded bracket solve uses
    /// the exact caloric source gain, never a numerical tangent or second heat
    /// march. Accepted external heat uses the original time-mean temperature,
    /// and the actual pressure/remap independently checks Delta H+W=Q_external.
    ///
    /// This is the existing uniform-temperature, pure-liquid source split,
    /// not a spatial heat solve or general viscous/free-surface integrator.
    /// The caller must justify uniform temperature and fixed total G over the
    /// interval (for example an unchanged planar heated area). Changing side
    /// area, reservoir capacity, radiation and changing source coefficients
    /// require their own coupled owners. A smaller timestep cannot establish
    /// the physical validity of spatial lumping. Repeated accepted steps retain
    /// the native affine motion and exact mass/enthalpy/momentum remap.
    ///
    /// The original chart, Cp, density, pressure/reference convention and
    /// endpoint transport sources stay intact. Fixed-pressure h already
    /// includes ambient p*dV; it is not charged again. Numerical projection
    /// loss and unresolved kinetic mapping changes never become heat.
    pub fn prepare_uniform_slab_bath_heat<'a>(
        &'a self,
        anchored_face: Face3,
        source: FreeSurfaceSlabBath,
        controls: FreeSurfaceSlabHeatControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabBathHeatPlan<'a, 's, 'm>, FreeSurfaceSlabBathHeatError> {
        checkpoint(gate, cx)?;
        if !source.reservoir_temperature_k.is_finite()
            || source.reservoir_temperature_k <= 0.0
            || !source.conductance_w_k.is_finite()
            || source.conductance_w_k < 0.0
            || !source.applied_power_w.is_finite()
        {
            return Err(invalid(
                "finite positive bath temperature and nonnegative conductance required",
            ));
        }
        let prepared = self.prepare_uniform_slab_heat_state(anchored_face, controls, gate, cx)?;
        let dt = Dd::from_f64(prepared.dt);
        let applied_heat = Dd::from_f64(source.applied_power_w) * dt;
        let conductance_time = Dd::from_f64(source.conductance_w_k) * dt;
        let temperature_drive = Dd::from_f64(source.reservoir_temperature_k)
            - Dd::from_f64(prepared.initial.temperature_k());
        let initial_drive = applied_heat + conductance_time * temperature_drive;
        let exposure = value(conductance_time / prepared.mass);
        if !value(applied_heat).is_finite()
            || !value(initial_drive).is_finite()
            || !exposure.is_finite()
            || (source.conductance_w_k > 0.0 && exposure == 0.0)
            || (source.applied_power_w != 0.0 && value(applied_heat) == 0.0)
        {
            return Err(invalid(
                "finite bath exposure or applied heat is unrepresentable",
            ));
        }
        let mut accepted_exchange = None;
        let solution = balance::solve_residual(
            prepared.mass,
            value(initial_drive / prepared.mass),
            controls.net_specific_enthalpy_bracket_j_kg,
            controls.maximum_heat_defect_j,
            controls.maximum_source_evaluations,
            |increment| -> Result<balance::EnergyResidual, FreeSurfaceSlabBathHeatError> {
                let mechanical = prepared.pressure_response(increment, gate, cx)?;
                let response = self.caloric.driven_heat_exchange_response_with_checkpoint(
                    prepared.initial,
                    value((initial_drive - mechanical.work_j) / prepared.mass),
                    exposure,
                    || gate.is_requested() || cx.is_cancel_requested(),
                )?;
                if response.exchange.state.phase() != SolidLiquidPhase::Liquid {
                    return Err(invalid(
                        "finite bath exchange leaves the original liquid phase",
                    ));
                }
                let residual = prepared.mass
                    * (Dd::from_f64(increment)
                        - Dd::from_f64(response.exchange.specific_heat_j_kg));
                let derivative = prepared.mass
                    + Dd::from_f64(response.jacobian.source_heat_gain)
                        * mechanical.work_derivative_kg;
                accepted_exchange = Some((increment.to_bits(), response.exchange));
                Ok(balance::EnergyResidual {
                    residual_j: residual,
                    derivative_kg: derivative,
                    work_j: mechanical.work_j,
                })
            },
        )
        .map_err(|failure| match failure {
            balance::Failure::Invalid(what) => invalid(what),
            balance::Failure::Evaluation(error) => error,
            balance::Failure::Balance {
                evaluations,
                residual_j,
            } => FreeSurfaceSlabHeatError::Balance {
                evaluations,
                residual_j,
            }
            .into(),
        })?;
        let (_, exchange) = accepted_exchange
            .filter(|(bits, _)| *bits == solution.increment_j_kg.to_bits())
            .ok_or_else(|| invalid("finite bath root lost its original caloric exchange"))?;
        // Relative time-mean T retains weak transfer near a large temperature datum.
        let reservoir_heat = conductance_time
            * (temperature_drive - Dd::from_f64(exchange.mean_temperature_change_k));
        let external = applied_heat + reservoir_heat;
        let defect = value(
            prepared.mass * Dd::from_f64(solution.increment_j_kg) + Dd::from_f64(solution.work_j)
                - external,
        );
        if !value(external).is_finite()
            || !value(reservoir_heat).is_finite()
            || !defect.is_finite()
            || defect.abs() > controls.maximum_heat_defect_j
        {
            return Err(invalid(
                "finite bath source and native pressure-work root do not close",
            ));
        }
        let duration_s = prepared.dt;
        let native = prepared.finish(solution, value(external), gate, cx)?;
        Ok(FreeSurfaceSlabBathHeatPlan {
            native,
            source,
            exchange,
            duration_s,
            applied_heat,
            reservoir_heat,
        })
    }
}

fn invalid(what: &'static str) -> FreeSurfaceSlabBathHeatError {
    FreeSurfaceSlabHeatError::Invalid(what).into()
}
