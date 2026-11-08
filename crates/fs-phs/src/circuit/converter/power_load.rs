//! Regulated electrical demand on the original converter bus coordinates.
//!
//! Above an authored foldback voltage, a load draws I=P/V until its current
//! capability binds. Below foldback it follows an explicit ohmic law, avoiding
//! the singular constant-power demand at zero volts. The resulting negative
//! incremental conductance in constant-power operation participates in the
//! SAME implicit converter solve; no previous-voltage current command is used.
//!
//! Load work leaves the converter through an external power boundary. It is
//! never deposited as converter conduction heat. Original L/C and all five
//! internal loss owners remain unchanged. This is a quasi-static regulated
//! input law, not a controller-bandwidth, efficiency, startup hysteresis or
//! semiconductor model. Calibration and voltage-domain applicability belong
//! to the caller. Both voltage polarities draw nonnegative load power.
//!
//! Evolution uses the original analytic discrete-gradient/Newton owner. Work
//! and heat use its discrete effort and midpoint quadrature, not an exact
//! continuous RMS or voltage-time integral. Nonlinear roots need not be unique;
//! convergence and energy closure do not certify stability or temporal error.

use super::{AveragedBuckBoost, ConductionHeat, ConverterStep, invalid, poll, positive};
use crate::{
    NonlinearPortPower, PhsError, PortHamiltonian, PreparedStepError, SolveBudget, StepWorkspace,
};
use fs_math::{dd::Dd, det};

/// An explicit continuous, bipolar, current-limited electrical input law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegulatedPowerLoad {
    requested_power_w: f64,
    maximum_current_a: f64,
    foldback_voltage_v: f64,
    maximum_voltage_v: f64,
    foldback_conductance_s: f64,
}

/// Active branch of the authored input law at one discrete electrical effort.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PowerLoadMode {
    Disabled,
    VoltageFoldback,
    CurrentLimited,
    ConstantPower,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerLoadPoint {
    pub mode: PowerLoadMode,
    pub voltage_v: f64,
    /// Signed current leaving the converter bus [A].
    pub current_a: f64,
    pub requested_power_w: f64,
    /// Actual V*I at this same effort, nonnegative for both polarities [W].
    pub delivered_power_w: f64,
    /// Nonnegative requested-minus-delivered power; multiplication roundoff at
    /// exact regulation can make actual power exceed the request by a few ulps.
    pub curtailed_power_w: f64,
    /// Branch derivative dI/dV, negative in constant-power operation [S].
    /// At a current-limit equality the capped derivative is selected; at the
    /// voltage-foldback equality the foldback derivative is selected.
    pub incremental_conductance_s: f64,
}
impl PowerLoadPoint {
    fn disabled(voltage_v: f64) -> Self {
        Self {
            mode: PowerLoadMode::Disabled,
            voltage_v,
            current_a: 0.0,
            requested_power_w: 0.0,
            delivered_power_w: 0.0,
            curtailed_power_w: 0.0,
            incremental_conductance_s: 0.0,
        }
    }
}

impl RegulatedPowerLoad {
    /// For |V|<=V_fold, the unconstrained law is I=P*V/V_fold²; above it,
    /// I=P/V. Its magnitude is then limited to I_max. Both transitions are
    /// continuous. Zero requested power disables demand without changing storage.
    /// No cutoff hysteresis or discontinuous turn-on is inferred.
    pub fn try_new(
        requested_power_w: f64,
        maximum_current_a: f64,
        foldback_voltage_v: f64,
        maximum_voltage_v: f64,
    ) -> Result<Self, PhsError> {
        let conductance = (requested_power_w / foldback_voltage_v) / foldback_voltage_v;
        if !requested_power_w.is_finite()
            || requested_power_w < 0.0
            || !positive(maximum_current_a)
            || !positive(foldback_voltage_v)
            || !positive(maximum_voltage_v)
            || maximum_voltage_v < foldback_voltage_v
            || !conductance.is_finite()
            || (requested_power_w > 0.0 && conductance <= 0.0)
        {
            return Err(invalid(
                "power load requires finite nonnegative demand and positive representable current/voltage limits",
            ));
        }
        Ok(Self {
            requested_power_w,
            maximum_current_a,
            foldback_voltage_v,
            maximum_voltage_v,
            foldback_conductance_s: conductance,
        })
    }
    pub const fn requested_power_w(&self) -> f64 {
        self.requested_power_w
    }
    pub const fn maximum_current_a(&self) -> f64 {
        self.maximum_current_a
    }
    pub const fn foldback_voltage_v(&self) -> f64 {
        self.foldback_voltage_v
    }
    pub const fn maximum_voltage_v(&self) -> f64 {
        self.maximum_voltage_v
    }

    /// Evaluate actual current, delivered power and the exact active-branch
    /// tangent. Nonfinite voltage, unsupported voltage, or an unrepresentable
    /// nonzero current/power refuses; demand is never extrapolated or fabricated.
    pub fn evaluate(&self, voltage_v: f64) -> Result<PowerLoadPoint, PhsError> {
        let magnitude = voltage_v.abs();
        if !voltage_v.is_finite() || magnitude > self.maximum_voltage_v {
            return Err(invalid(
                "regulated power load voltage is outside its authored domain",
            ));
        }
        if self.requested_power_w == 0.0 {
            return Ok(PowerLoadPoint::disabled(voltage_v));
        }
        let (unlimited, slope, branch) = if magnitude <= self.foldback_voltage_v {
            (
                (self.requested_power_w / self.foldback_voltage_v)
                    * (magnitude / self.foldback_voltage_v),
                self.foldback_conductance_s,
                PowerLoadMode::VoltageFoldback,
            )
        } else {
            let current = self.requested_power_w / magnitude;
            (
                current,
                -(current / magnitude),
                PowerLoadMode::ConstantPower,
            )
        };
        let (current, derivative, mode) = if unlimited >= self.maximum_current_a {
            (self.maximum_current_a, 0.0, PowerLoadMode::CurrentLimited)
        } else {
            (unlimited, slope, branch)
        };
        let current_a = current.copysign(voltage_v);
        let delivered = voltage_v * current_a;
        if !current_a.is_finite()
            || !delivered.is_finite()
            || !derivative.is_finite()
            || (magnitude > 0.0 && (current == 0.0 || delivered == 0.0))
            || (delivered > self.requested_power_w
                && (delivered - self.requested_power_w) / self.requested_power_w
                    > 16.0 * f64::EPSILON)
        {
            return Err(invalid(
                "regulated load current, work or tangent is unrepresentable",
            ));
        }
        Ok(PowerLoadPoint {
            mode,
            voltage_v,
            current_a,
            requested_power_w: self.requested_power_w,
            delivered_power_w: delivered,
            curtailed_power_w: (self.requested_power_w - delivered).max(0.0),
            incremental_conductance_s: derivative,
        })
    }
}

/// Explicit original Newton/backtracking work and independent energy allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerLoadControls {
    pub solve: SolveBudget,
    /// At most 16 halvings per Newton update; zero disables backtracking.
    pub backtracking_halvings: u8,
    pub absolute_energy_tolerance_j: f64,
    /// Relative to this interval's actual work/heat/storage change, never a
    /// large untouched capacitor-energy datum. Must lie in [0,1).
    pub relative_energy_tolerance: f64,
}

#[derive(Clone, Debug)]
pub struct PowerLoadedStep {
    /// Native converter record; supply is NET source work minus exported load
    /// work. Dissipation and its five channels contain converter losses only.
    pub converter: ConverterStep,
    pub loads: [PowerLoadPoint; 2],
    pub source_work_j: [f64; 2],
    /// Delivered load energy, exported once through the electrical boundary.
    pub load_work_j: [f64; 2],
    /// Numerical midpoint curtailment, not the continuous controller's energy.
    pub curtailed_load_energy_j: [f64; 2],
}

/// Original three-state converter plus one reusable nonlinear solve workspace.
/// Loads are held commands passed to each step, so changing demand never resets
/// the converter's original flux/charges or introduces parameter-storage work.
pub struct PreparedPowerLoadedBuckBoost {
    model: AveragedBuckBoost,
    workspace: StepWorkspace,
    inverse_storage: [f64; 3],
    duration_s: f64,
    controls: PowerLoadControls,
}

impl AveragedBuckBoost {
    /// Prepare the existing three-state analytic implicit owner. Allocation is
    /// cold; successful warm steps and both physical callbacks allocate nothing.
    /// Nonlinear Jacobian factorization remains bounded by the supplied Newton
    /// budget. This is a numerical averaged converter with quasi-static loads.
    pub fn prepare_power_loaded_controlled(
        self,
        duration_s: f64,
        controls: PowerLoadControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPowerLoadedBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid(
                "power-loaded converter requires positive duration and finite energy allowances",
            )
            .into());
        }
        controls.solve.validate()?;
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        let inverse_storage = [
            1.0 / self.parameters.inductance_h,
            1.0 / self.parameters.bus_capacitance_f[0],
            1.0 / self.parameters.bus_capacitance_f[1],
        ];
        poll(&mut cancelled)?;
        Ok(PreparedPowerLoadedBuckBoost {
            model: self,
            workspace,
            inverse_storage,
            duration_s,
            controls,
        })
    }
}

pub(super) fn load_point(
    load: Option<RegulatedPowerLoad>,
    voltage: f64,
) -> Result<PowerLoadPoint, PhsError> {
    if !voltage.is_finite() {
        return Err(invalid("power-load effort is not finite"));
    }
    load.map_or_else(
        || Ok(PowerLoadPoint::disabled(voltage)),
        |load| load.evaluate(voltage),
    )
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}

impl PreparedPowerLoadedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Solve both load currents and the full converter together. External
    /// terminal currents are held and positive INTO each bus. A missing load
    /// exports zero work; positive demand draws power at either voltage polarity.
    /// Output voltages are the native DISCRETE efforts, not exact time averages.
    /// Domain checks cover initial/final and discrete voltages; they do not
    /// certify unseen continuous-time overshoots. Failure preserves both outputs.
    #[allow(clippy::too_many_arguments)]
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; 3],
        discrete_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        for p in 0..2 {
            load_point(loads[p], initial[p + 1] * self.inverse_storage[p + 1])?;
        }
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if effort.len() != 3 || out.len() != 3 {
                return None;
            }
            out.fill(0.0);
            let mut power = [0.0; 2];
            for p in 0..2 {
                let point = load_point(loads[p], effort[p + 1]).ok()?;
                out[p + 1] = point.current_a;
                power[p] = point.delivered_power_w;
            }
            // The load is an external work receiver, not converter heat.
            Some(NonlinearPortPower {
                supplied_power_w: -sum(power),
                dissipated_power_w: 0.0,
            })
        };
        let inverse_storage = self.inverse_storage;
        let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
            if direction.len() != 3 || out.len() != 3 {
                return false;
            }
            for i in 0..3 {
                out[i] = inverse_storage[i] * direction[i];
            }
            out.iter().all(|v| v.is_finite())
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if effort.len() != 3 || de.len() != 3 || out.len() != 3 {
                return false;
            }
            out.fill(0.0);
            for p in 0..2 {
                let Ok(point) = load_point(loads[p], effort[p + 1]) else {
                    return false;
                };
                out[p + 1] = point.incremental_conductance_s * de[p + 1];
            }
            out.iter().all(|v| v.is_finite())
        };
        let (mut candidate, mut voltage) = ([0.0; 3], [0.0; 2]);
        let native = self.workspace.step_into_supplied_controlled(
            &self.model.system,
            initial,
            terminal_current_a,
            self.duration_s,
            &mut candidate,
            &mut voltage,
            &port,
            Some(&hessian),
            Some(&tangent),
            &mut cancelled,
        )?;
        let effort = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("power load requires the accepted native discrete effort"))?;
        let mut points = [PowerLoadPoint::disabled(0.0); 2];
        for p in 0..2 {
            load_point(loads[p], candidate[p + 1] * self.inverse_storage[p + 1])?;
            points[p] = load_point(loads[p], effort[p + 1])?;
        }
        let load_work_j: [f64; 2] = points.map(|p| self.duration_s * p.delivered_power_w);
        let curtailed_load_energy_j: [f64; 2] =
            points.map(|p| self.duration_s * p.curtailed_power_w);
        let source_work_j: [f64; 2] =
            core::array::from_fn(|p| (self.duration_s * voltage[p]) * terminal_current_a[p]);
        let terminal_work_j: [f64; 2] =
            core::array::from_fn(|p| sum([source_work_j[p], -load_work_j[p]]));
        let joules = self.model.effective_resistance_ohm.map(|r| {
            let integrated_amplitude = (det::sqrt(r) * effort[0]) * det::sqrt(self.duration_s);
            integrated_amplitude * integrated_amplitude
        });
        let scale = sum(source_work_j.map(f64::abs))
            .max(sum(load_work_j))
            .max(native.delta_h.abs())
            .max(native.dissipated);
        let allowance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            sum([native.delta_h, native.dissipated, -sum(terminal_work_j)]),
            sum([native.supplied, -sum(terminal_work_j)]),
            sum(joules.into_iter().chain([-native.dissipated])),
        ];
        if !scale.is_finite()
            || !allowance.is_finite()
            || candidate
                .iter()
                .chain(&voltage)
                .chain(&joules)
                .chain(&load_work_j)
                .chain(&source_work_j)
                .chain(&curtailed_load_energy_j)
                .chain(&terminal_work_j)
                .any(|v| !v.is_finite())
            || residuals
                .iter()
                .any(|r| !r.is_finite() || r.abs() > allowance)
        {
            return Err(invalid(
                "power-loaded converter source, exported work and conduction heat do not close",
            )
            .into());
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        discrete_bus_voltage_v.copy_from_slice(&voltage);
        Ok(PowerLoadedStep {
            converter: ConverterStep {
                native,
                terminal_work_j,
                conduction: ConductionHeat { joules },
            },
            loads: points,
            source_work_j,
            load_work_j,
            curtailed_load_energy_j,
        })
    }
}
