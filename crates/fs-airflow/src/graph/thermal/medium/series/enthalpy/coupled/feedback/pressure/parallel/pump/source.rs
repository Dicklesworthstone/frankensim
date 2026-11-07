//! Electrical source limits on calibrated, material-dependent pump operation.
//!
//! A speed trial prepares the original calibrated head AND independent input
//! power curves, then solves the complete hydraulic/material operating point.
//! Electrical power is observed at that actual flow; no guessed-flow cutback,
//! affinity scaling, efficiency fit or clipping of reported power is applied.
//!
//! A feasible command is retained. Otherwise the supplied lower speed must be
//! feasible and brackets the command's excess draw. Bisection retains that
//! feasible/infeasible bracket and accepts only an actual pressure-closed point
//! below the source cap and within its stated power tolerance. It claims no
//! globally maximal speed, monotone response, uniqueness or stability. A branch
//! jump, unsupported trial or unresolved power boundary refuses atomically.
//!
//! This is held quasi-steady speed selection. It adds no rotor acceleration
//! energy, source inventory, controller dynamics or finite-reservoir timeline.
//! The selected curve keeps its original loss destinations; hydraulic work is
//! never relabeled as pump heat. Changed walls/inlet/history need a new solve.

use super::super::solid::CoupledParallelLiquidSolution;
use super::solid::PumpCouplingControls;
use super::speed::PumpSpeedMap;
use super::{
    LiquidFeedbackBudget, LiquidResistanceModel, ParallelLiquidNetwork, ParallelLiquidSolution,
    ParallelPressureBudget, ParallelPressureError, PumpHeadCurve, PumpOperatingControls,
    PumpOperatingError, PumpOperatingPoint, PumpSupplyState, checkpoint,
};
use crate::AirflowError;
use crate::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;

/// Cold speed trials are bounded independently of each native pressure solve.
pub const MAX_POWER_SPEED_EVALUATIONS: usize = 64;

/// Source capability and numerical resolution of a bracketed power boundary.
#[derive(Clone, Copy, Debug)]
pub struct PumpPowerLimitControls {
    /// Calibrated lower speed [rev/min], no greater than the requested command.
    /// It must have an admitted flowing operating point if cutback is needed.
    pub lower_speed_rpm: f64,
    /// Nonnegative electrical draw capability [W], never a hydraulic-power cap.
    pub maximum_electrical_input_w: f64,
    /// Maximum shortfall below the cap accepted after curtailment [W].
    /// The actual draw must still be <= the cap, without an overshoot allowance.
    pub absolute_power_tolerance_w: f64,
    /// Additional relative shortfall allowance in [0,1), applied to the cap.
    pub relative_power_tolerance: f64,
    /// Full fixed-speed pump operating points, including command and lower end.
    /// Must lie in 2..=MAX_POWER_SPEED_EVALUATIONS.
    pub max_speed_evaluations: usize,
    /// Total original network pressure candidates across ALL speed trials.
    /// Must lie in 2..=8192; each native pressure-solve cap also remains active.
    pub max_total_network_evaluations: usize,
}

impl PumpPowerLimitControls {
    fn admit(self, map: &PumpSpeedMap, command: f64) -> Result<(), PumpPowerLimitError> {
        let (lower, upper) = map.speed_domain_rpm();
        if !self.lower_speed_rpm.is_finite()
            || !command.is_finite()
            || self.lower_speed_rpm < lower
            || command > upper
            || self.lower_speed_rpm > command
            || !self.maximum_electrical_input_w.is_finite()
            || self.maximum_electrical_input_w < 0.0
            || !self.absolute_power_tolerance_w.is_finite()
            || self.absolute_power_tolerance_w < 0.0
            || !(0.0..1.0).contains(&self.relative_power_tolerance)
            || !(2..=MAX_POWER_SPEED_EVALUATIONS).contains(&self.max_speed_evaluations)
            || !(2..=8192).contains(&self.max_total_network_evaluations)
        {
            return Err(PumpPowerLimitError::Invalid(
                "power control needs an admitted speed interval, finite source/tolerances and complete work caps",
            ));
        }
        Ok(())
    }

    fn shortfall_tolerance(self) -> f64 {
        self.absolute_power_tolerance_w
            .max(self.relative_power_tolerance * self.maximum_electrical_input_w)
    }
}

/// Owned calibrated curve and the original accepted state at that same speed.
/// No self-reference or detached approximation to the hydraulic state is stored.
#[derive(Debug)]
pub struct PowerLimitedPumpOperatingPoint<T> {
    pub curve: PumpHeadCurve,
    pub state: T,
    pub supply: PumpSupplyState,
    pub network_pressure_pa: f64,
    pub pressure_residual_pa: f64,
    pub commanded_speed_rpm: f64,
    pub maximum_electrical_input_w: f64,
    /// Nonnegative actual cap minus independent electrical consumption [W].
    pub electrical_power_margin_w: f64,
    pub curtailed: bool,
    pub speed_evaluations: usize,
    pub total_network_evaluations: usize,
    pub total_flow_evaluations: usize,
    pub total_solid_solves: usize,
}

/// No failed search publishes a partially matched network or solid state.
#[derive(Debug)]
pub enum PumpPowerLimitError {
    Invalid(&'static str),
    Operating(PumpOperatingError),
    NetworkBudget {
        performed: usize,
        maximum: usize,
    },
    /// The authored lower endpoint does not bracket the command's power excess.
    /// This is not a proof that all intermediate speeds are infeasible.
    NotBracketed {
        lower_speed_rpm: f64,
        lower_electrical_input_w: f64,
        commanded_electrical_input_w: f64,
        maximum_electrical_input_w: f64,
    },
    NotConverged {
        speed_evaluations: usize,
        feasible_speed_rpm: f64,
        infeasible_speed_rpm: f64,
        feasible_electrical_input_w: f64,
        maximum_electrical_input_w: f64,
    },
}

impl From<PumpOperatingError> for PumpPowerLimitError {
    fn from(error: PumpOperatingError) -> Self {
        Self::Operating(error)
    }
}
impl core::fmt::Display for PumpPowerLimitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "calibrated pump source limit: {self:?}")
    }
}
impl std::error::Error for PumpPowerLimitError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Operating(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Work {
    speeds: usize,
    networks: usize,
    flows: usize,
    solids: usize,
}

// Strip only the temporary curve borrow from the native result. The original
// actual network/solid state and all counters survive without cloning them.
struct Trial<T> {
    state: T,
    supply: PumpSupplyState,
    pressure_pa: f64,
    residual_pa: f64,
    networks: usize,
    flows: usize,
    solids: usize,
}
impl<T> From<PumpOperatingPoint<'_, T>> for Trial<T> {
    fn from(point: PumpOperatingPoint<'_, T>) -> Self {
        Self {
            state: point.state,
            supply: point.supply,
            pressure_pa: point.network_pressure_pa,
            residual_pa: point.pressure_residual_pa,
            networks: point.network_evaluations,
            flows: point.total_flow_evaluations,
            solids: point.total_solid_solves,
        }
    }
}

struct Candidate<T> {
    curve: PumpHeadCurve,
    trial: Trial<T>,
    speed_rpm: f64,
    electrical_input_w: f64,
}

fn check(cx: &Cx<'_>) -> Result<(), PumpPowerLimitError> {
    checkpoint(cx)
        .map_err(ParallelPressureError::from)
        .map_err(PumpOperatingError::from)?;
    Ok(())
}

fn add(a: usize, b: usize) -> Result<usize, PumpPowerLimitError> {
    a.checked_add(b).ok_or(PumpPowerLimitError::Invalid(
        "pump source work count overflow",
    ))
}

fn probe<T>(
    cx: &Cx<'_>,
    map: &PumpSpeedMap,
    speed: f64,
    source: PumpPowerLimitControls,
    used: &mut Work,
    evaluate: &mut impl FnMut(&PumpHeadCurve, Work) -> Result<Trial<T>, PumpPowerLimitError>,
) -> Result<Candidate<T>, PumpPowerLimitError> {
    check(cx)?;
    if used.speeds >= source.max_speed_evaluations {
        return Err(PumpPowerLimitError::Invalid(
            "pump speed probe exceeded its admitted count",
        ));
    }
    if source
        .max_total_network_evaluations
        .saturating_sub(used.networks)
        < 2
    {
        return Err(PumpPowerLimitError::NetworkBudget {
            performed: used.networks,
            maximum: source.max_total_network_evaluations,
        });
    }
    let curve = map.at_speed(speed)?;
    check(cx)?;
    let trial = evaluate(&curve, *used)?;
    let input = trial
        .supply
        .drive_power
        .ok_or(PumpPowerLimitError::Invalid(
            "source limiting requires the original independent electrical input curve",
        ))?
        .electrical_input_w;
    if !input.is_finite()
        || input < 0.0
        || trial.supply.speed_rpm != Some(speed)
        || curve.speed_rpm() != Some(speed)
        || trial.networks == 0
    {
        return Err(PumpPowerLimitError::Invalid(
            "source limit received an inconsistent calibrated pump operating point",
        ));
    }
    used.speeds = add(used.speeds, 1)?;
    used.networks = add(used.networks, trial.networks)?;
    used.flows = add(used.flows, trial.flows)?;
    used.solids = add(used.solids, trial.solids)?;
    if used.networks > source.max_total_network_evaluations {
        return Err(PumpPowerLimitError::NetworkBudget {
            performed: used.networks,
            maximum: source.max_total_network_evaluations,
        });
    }
    check(cx)?;
    Ok(Candidate {
        curve,
        trial,
        speed_rpm: speed,
        electrical_input_w: input,
    })
}

fn finish<T>(
    cx: &Cx<'_>,
    point: Candidate<T>,
    command: f64,
    source: PumpPowerLimitControls,
    used: Work,
) -> Result<PowerLimitedPumpOperatingPoint<T>, PumpPowerLimitError> {
    let margin = source.maximum_electrical_input_w - point.electrical_input_w;
    if !margin.is_finite() || margin < 0.0 {
        return Err(PumpPowerLimitError::Invalid(
            "accepted pump exceeds its electrical source cap",
        ));
    }
    check(cx)?;
    Ok(PowerLimitedPumpOperatingPoint {
        curve: point.curve,
        state: point.trial.state,
        supply: point.trial.supply,
        network_pressure_pa: point.trial.pressure_pa,
        pressure_residual_pa: point.trial.residual_pa,
        commanded_speed_rpm: command,
        maximum_electrical_input_w: source.maximum_electrical_input_w,
        electrical_power_margin_w: margin,
        curtailed: point.speed_rpm != command,
        speed_evaluations: used.speeds,
        total_network_evaluations: used.networks,
        total_flow_evaluations: used.flows,
        total_solid_solves: used.solids,
    })
}

fn match_power<T>(
    cx: &Cx<'_>,
    map: &PumpSpeedMap,
    command: f64,
    source: PumpPowerLimitControls,
    mut evaluate: impl FnMut(&PumpHeadCurve, Work) -> Result<Trial<T>, PumpPowerLimitError>,
) -> Result<PowerLimitedPumpOperatingPoint<T>, PumpPowerLimitError> {
    check(cx)?;
    source.admit(map, command)?;
    let mut used = Work::default();
    let commanded = probe(cx, map, command, source, &mut used, &mut evaluate)?;
    if commanded.electrical_input_w <= source.maximum_electrical_input_w {
        return finish(cx, commanded, command, source, used);
    }
    let commanded_power = commanded.electrical_input_w;
    drop(commanded);
    if source.lower_speed_rpm == command {
        return Err(PumpPowerLimitError::NotBracketed {
            lower_speed_rpm: command,
            lower_electrical_input_w: commanded_power,
            commanded_electrical_input_w: commanded_power,
            maximum_electrical_input_w: source.maximum_electrical_input_w,
        });
    }
    let mut lower = probe(
        cx,
        map,
        source.lower_speed_rpm,
        source,
        &mut used,
        &mut evaluate,
    )?;
    if lower.electrical_input_w > source.maximum_electrical_input_w {
        return Err(PumpPowerLimitError::NotBracketed {
            lower_speed_rpm: lower.speed_rpm,
            lower_electrical_input_w: lower.electrical_input_w,
            commanded_electrical_input_w: commanded_power,
            maximum_electrical_input_w: source.maximum_electrical_input_w,
        });
    }
    let mut high = command;
    loop {
        if source.maximum_electrical_input_w - lower.electrical_input_w
            <= source.shortfall_tolerance()
        {
            return finish(cx, lower, command, source, used);
        }
        if used.speeds == source.max_speed_evaluations {
            break;
        }
        let middle = 0.5 * lower.speed_rpm + 0.5 * high;
        if middle == lower.speed_rpm || middle == high {
            break;
        }
        let point = probe(cx, map, middle, source, &mut used, &mut evaluate)?;
        if point.electrical_input_w <= source.maximum_electrical_input_w {
            lower = point;
        } else {
            high = middle;
        }
    }
    Err(PumpPowerLimitError::NotConverged {
        speed_evaluations: used.speeds,
        feasible_speed_rpm: lower.speed_rpm,
        infeasible_speed_rpm: high,
        feasible_electrical_input_w: lower.electrical_input_w,
        maximum_electrical_input_w: source.maximum_electrical_input_w,
    })
}

fn remaining_flow(
    work: ParallelPressureBudget,
    used: Work,
) -> Result<ParallelPressureBudget, PumpPowerLimitError> {
    let remaining = work.max_flow_evaluations.saturating_sub(used.flows);
    if remaining < 2 {
        return Err(
            PumpOperatingError::Network(ParallelPressureError::FlowBudget {
                performed: used.flows,
                maximum: work.max_flow_evaluations,
            })
            .into(),
        );
    }
    Ok(ParallelPressureBudget {
        max_flow_evaluations: remaining,
        ..work
    })
}

fn remaining_network(
    operating: PumpOperatingControls,
    source: PumpPowerLimitControls,
    used: Work,
) -> PumpOperatingControls {
    PumpOperatingControls {
        max_network_evaluations: operating.max_network_evaluations.min(
            source
                .max_total_network_evaluations
                .saturating_sub(used.networks),
        ),
        ..operating
    }
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Source-limited speed at the actual fixed-wall, material-dependent network
    /// operating point. `work.max_flow_evaluations` spans ALL speed/pressure
    /// trials; source controls additionally bound speed and network evaluations.
    /// All native pressure/material/power/omitted-heating gates remain active.
    ///
    /// The returned owned curve is a held-speed calibration. Reusing it after
    /// changing the fluid/walls does not retain this operating point's source
    /// bound; run this solve again for the new state. No accepted history moves.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_power_limited_pump(
        &self,
        cx: &Cx<'_>,
        walls: &[f64],
        map: &PumpSpeedMap,
        commanded_speed_rpm: f64,
        source: PumpPowerLimitControls,
        local: LiquidFeedbackBudget,
        work: ParallelPressureBudget,
        operating: PumpOperatingControls,
    ) -> Result<PowerLimitedPumpOperatingPoint<ParallelLiquidSolution<'a>>, PumpPowerLimitError>
    {
        check(cx)?;
        operating.admit()?;
        work.admit().map_err(PumpOperatingError::from)?;
        match_power(cx, map, commanded_speed_rpm, source, |curve, used| {
            let point = self.solve_with_pump(
                cx,
                walls,
                curve,
                local,
                remaining_flow(work, used)?,
                remaining_network(operating, source, used),
            )?;
            Ok(point.into())
        })
    }

    /// The same source limiter around original common-solid/material/branch
    /// coupling. Every speed and pressure trial re-solves the complete solid.
    /// Existing total flow and solid caps in `controls` span the ENTIRE speed
    /// search, and every original inner convergence gate remains authoritative.
    /// Callbacks must be disposable evaluations from the same accepted history;
    /// no failed or superseded candidate may publish external state.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_power_limited_pump_coupled_iqn<R, F>(
        &self,
        cx: &Cx<'_>,
        initial_walls_k: &[f64],
        map: &PumpSpeedMap,
        commanded_speed_rpm: f64,
        source: PumpPowerLimitControls,
        controls: &PumpCouplingControls,
        mut solid: F,
    ) -> Result<
        PowerLimitedPumpOperatingPoint<CoupledParallelLiquidSolution<'a, R>>,
        PumpPowerLimitError,
    >
    where
        F: FnMut(
            &Cx<'_>,
            &[CorrelationRobinBoundary],
        ) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        check(cx)?;
        controls.operating_point.admit()?;
        controls
            .network
            .work
            .admit()
            .map_err(PumpOperatingError::from)?;
        if controls.max_total_solid_solves == 0 {
            return Err(PumpPowerLimitError::Invalid(
                "whole-search solid callback budget must be positive",
            ));
        }
        match_power(cx, map, commanded_speed_rpm, source, |curve, used| {
            let mut inner = *controls;
            inner.network.work = remaining_flow(controls.network.work, used)?;
            inner.operating_point = remaining_network(controls.operating_point, source, used);
            inner.max_total_solid_solves =
                controls.max_total_solid_solves.saturating_sub(used.solids);
            if inner.max_total_solid_solves == 0 {
                return Err(PumpOperatingError::SolidBudget {
                    performed: used.solids,
                    maximum: controls.max_total_solid_solves,
                }
                .into());
            }
            let point =
                self.solve_with_pump_coupled_iqn(cx, initial_walls_k, curve, &inner, &mut solid)?;
            Ok(point.into())
        })
    }
}
