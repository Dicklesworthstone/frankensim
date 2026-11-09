//! Reversible nonlinear winding storage in the original averaged power stage.
//!
//! The state remains `[flux linkage, input charge, output charge]`. The complete
//! source integral `H = integral I dλ` replaces the linear magnetic potential;
//! the two original bridge exchanges, finite buses and five resistance owners
//! are unchanged. Saturation affects current, magnetic energy and differential
//! inductance together. This separate owner cannot enter a linear-only adapter
//! that would reconstruct the winding from a fixed L.
//!
//! Work and heat use the native discrete-gradient quadrature. Winding current
//! is the energy secant over flux change, not an exact time-average or RMS value;
//! bus efforts are midpoint voltages. Conduction heat `dt R I_discrete²` excludes
//! unresolved PWM ripple. Source interpolation and time-step accuracy remain
//! separate from exact constitutive work. No hysteresis, switching, core loss,
//! temperature law or source extrapolation is added by this reversible model.

use super::{
    AveragedBuckBoost, BridgeModulation, CONDUCTION_CHANNELS, ConductionHeat, invalid, poll,
};
use crate::circuit::{
    CircuitError, capacitance::storage_with_inductor_curve, inductance::InductorFluxCurve,
};
use crate::{
    PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage, SolveBudget,
    StepWorkspace,
};
use fs_math::dd::Dd;

/// Regulated demand solved together with saturated winding and both finite buses.
pub mod power_load;
pub use power_load::NonlinearPowerLoadedStep;

/// Explicit native solve limits and independent physical-energy allowances.
/// Defaults are strict: no allowance is raised to accept a trial. Small charge
/// increments on a charged bus can require an explicit absolute allowance for
/// endpoint representation (one half charge ULP times its discrete voltage),
/// independently of the source interpolation and time-discretization error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearConverterControls {
    pub solve: SolveBudget,
    /// Existing native globalization, 0..=16; unsupported source probes refuse.
    pub backtracking_halvings: u8,
    pub absolute_energy_tolerance_j: f64,
    /// Applied to actual transfers, never to unchanged background storage [0,1).
    pub relative_balance_tolerance: f64,
}

impl Default for NonlinearConverterControls {
    fn default() -> Self {
        Self {
            solve: SolveBudget::default(),
            backtracking_halvings: 0,
            absolute_energy_tolerance_j: 0.0,
            relative_balance_tolerance: 1e-12,
        }
    }
}

/// The original three-state, two-port stage with a bounded winding potential.
pub struct NonlinearAveragedBuckBoost {
    bus_capacitance_f: [f64; 2],
    modulation: BridgeModulation,
    effective_resistance_ohm: [f64; CONDUCTION_CHANNELS],
    curve: InductorFluxCurve,
    system: PortHamiltonian,
}

/// Independently evaluated physical accounts on one native discrete step.
/// `conduction` partitions `native.dissipated`; they are not separate debits.
#[derive(Clone, Debug)]
pub struct NonlinearConverterStep {
    pub native: PreparedStepRecord,
    pub magnetic_energy_change_j: f64,
    /// Input and output capacitor energy changes, respectively [J].
    pub bus_energy_change_j: [f64; 2],
    /// Signed work INTO each original bus: `dt I_terminal V_discrete` [J].
    pub terminal_work_j: [f64; 2],
    pub conduction: ConductionHeat,
    /// Exact magnetic-energy secant, also used in both bridge current laws [A].
    pub discrete_current_a: f64,
}

impl AveragedBuckBoost {
    /// Replace only the original winding potential, retaining the SAME J/R/G.
    /// The construction-time linear L is discarded, with no duplicate energy.
    /// Use the returned owner's `state` to initialize from physical current;
    /// changing an already energized winding law requires external parameter work.
    /// The complete curve survives `into_system` and native interconnection.
    pub fn with_inductor_curve(
        self,
        curve: InductorFluxCurve,
    ) -> Result<NonlinearAveragedBuckBoost, CircuitError> {
        let bus_capacitance_f = self.parameters.bus_capacitance_f;
        let mut system = self.system;
        system.storage = storage_with_inductor_curve(
            QuadraticStorage::new(
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0 / bus_capacitance_f[0],
                    0.0,
                    0.0,
                    0.0,
                    1.0 / bus_capacitance_f[1],
                ],
                3,
            )?,
            0,
            &curve,
        )?;
        Ok(NonlinearAveragedBuckBoost {
            bus_capacitance_f,
            modulation: self.modulation,
            effective_resistance_ohm: self.effective_resistance_ohm,
            curve,
            system,
        })
    }
}

impl NonlinearAveragedBuckBoost {
    pub const fn bus_capacitance_f(&self) -> [f64; 2] {
        self.bus_capacitance_f
    }
    pub const fn modulation(&self) -> BridgeModulation {
        self.modulation
    }
    pub const fn effective_resistance_ohm(&self) -> &[f64; CONDUCTION_CHANNELS] {
        &self.effective_resistance_ohm
    }
    pub const fn flux_curve(&self) -> &InductorFluxCurve {
        &self.curve
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    /// Transfer the complete nonlinear owner into native composition.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Original flux and charges at a source-supported signed winding current.
    pub fn state(&self, current_a: f64, bus_voltage_v: [f64; 2]) -> Result<[f64; 3], CircuitError> {
        let state = [
            self.curve.flux_wb(current_a)?,
            self.bus_capacitance_f[0] * bus_voltage_v[0],
            self.bus_capacitance_f[1] * bus_voltage_v[1],
        ];
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid("nonlinear converter state or energy is unrepresentable").into());
        }
        Ok(state)
    }

    /// Prepare the native nonlinear solve with exact source-owned derivatives.
    /// Preparation allocates fixed three-state scratch; warm steps and source
    /// callbacks allocate nothing. Duty and terminal currents are held per step.
    pub fn prepare_controlled(
        self,
        duration_s: f64,
        controls: NonlinearConverterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedNonlinearAveragedBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_balance_tolerance)
        {
            return Err(invalid(
                "nonlinear converter needs positive duration, nonnegative absolute allowance and relative allowance in [0,1)",
            )
            .into());
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedNonlinearAveragedBuckBoost {
            converter: self,
            workspace,
            duration_s,
            controls,
        })
    }
}

/// Fixed-duty nonlinear propagator; all accepted state belongs to the caller.
pub struct PreparedNonlinearAveragedBuckBoost {
    converter: NonlinearAveragedBuckBoost,
    workspace: StepWorkspace,
    duration_s: f64,
    controls: NonlinearConverterControls,
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}
fn product(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::from_f64(1.0), |p, v| p * Dd::from_f64(v))
        .to_f64()
}

impl PreparedNonlinearAveragedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        self.converter.system()
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn flux_curve(&self) -> &InductorFluxCurve {
        self.converter.flux_curve()
    }

    /// Advance original flux and both charges, using the work-conjugate winding
    /// current in the reciprocal voltage/current equations. Output voltages are
    /// native discrete efforts, not exact continuous-time interval averages.
    /// Every source, solve, account or cancellation refusal preserves BOTH
    /// output buffers; the workspace has no accepted history and is retryable.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[f64; 2],
        next: &mut [f64; 3],
        discrete_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearConverterStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let (mut candidate, mut voltage) = ([0.0; 3], [0.0; 2]);
        let native = self.workspace.step_into_controlled(
            &self.converter.system,
            initial,
            terminal_current_a,
            self.duration_s,
            &mut candidate,
            &mut voltage,
            &mut cancelled,
        )?;
        let terminal_work_j: [f64; 2] =
            core::array::from_fn(|i| product([self.duration_s, terminal_current_a[i], voltage[i]]));
        let record = self.account_interval(
            initial,
            &candidate,
            &voltage,
            terminal_work_j,
            terminal_work_j.map(f64::abs),
            native,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        discrete_bus_voltage_v.copy_from_slice(&voltage);
        Ok(record)
    }

    // Both ordinary and regulated-load steps retain these SAME original source,
    // capacitor, bridge and five conduction accounts. Terminal scales include
    // actual source/export work before cancellation, never background storage.
    fn account_interval(
        &self,
        initial: &[f64; 3],
        candidate: &[f64; 3],
        voltage: &[f64; 2],
        terminal_work_j: [f64; 2],
        terminal_work_scale_j: [f64; 2],
        native: PreparedStepRecord,
    ) -> Result<NonlinearConverterStep, PreparedStepError> {
        let c = &self.converter;
        let h = self.duration_s;
        let magnetic_energy_change_j = c
            .curve
            .energy_change_j(initial[0], candidate[0])
            .map_err(|_| invalid("nonlinear converter magnetic work is outside source support"))?;
        let discrete_current_a = self
            .workspace
            .last_step_discrete_effort()
            .and_then(|effort| effort.first())
            .copied()
            .ok_or_else(|| invalid("nonlinear converter has no accepted winding effort"))?;
        let bus_energy_change_j: [f64; 2] = core::array::from_fn(|i| {
            let before = Dd::from_f64(initial[i + 1]);
            let after = Dd::from_f64(candidate[i + 1]);
            ((after - before) * (Dd::from_f64(0.5) * before + Dd::from_f64(0.5) * after)
                / Dd::from_f64(c.bus_capacitance_f[i]))
            .to_f64()
        });
        let joules = c
            .effective_resistance_ohm
            .map(|r| product([h, r, discrete_current_a, discrete_current_a]));
        // Signed internal work INTO the winding. The same quantities leave the
        // two buses; neither ideal bridge consumes heat or supplies net energy.
        let exchange = [
            product([
                h,
                c.modulation.input_high_fraction,
                voltage[0],
                discrete_current_a,
            ]),
            product([
                h,
                -c.modulation.output_high_fraction,
                voltage[1],
                discrete_current_a,
            ]),
        ];
        let heat = sum(joules);
        let supplied = sum(terminal_work_j);
        let stored = sum([
            magnetic_energy_change_j,
            bus_energy_change_j[0],
            bus_energy_change_j[1],
        ]);
        let stored_scale = sum([
            magnetic_energy_change_j.abs(),
            bus_energy_change_j[0].abs(),
            bus_energy_change_j[1].abs(),
        ]);
        let supply_scale = sum(terminal_work_scale_j);
        let accounts = [
            (
                sum([magnetic_energy_change_j, heat, -exchange[0], -exchange[1]]),
                sum([
                    magnetic_energy_change_j.abs(),
                    heat,
                    exchange[0].abs(),
                    exchange[1].abs(),
                ]),
            ),
            (
                sum([bus_energy_change_j[0], exchange[0], -terminal_work_j[0]]),
                sum([
                    bus_energy_change_j[0].abs(),
                    exchange[0].abs(),
                    terminal_work_scale_j[0],
                ]),
            ),
            (
                sum([bus_energy_change_j[1], exchange[1], -terminal_work_j[1]]),
                sum([
                    bus_energy_change_j[1].abs(),
                    exchange[1].abs(),
                    terminal_work_scale_j[1],
                ]),
            ),
            (sum([stored, -native.delta_h]), stored_scale),
            (sum([heat, -native.dissipated]), heat),
            (sum([supplied, -native.supplied]), supply_scale),
            (
                sum([stored, heat, -supplied]),
                sum([stored_scale, heat, supply_scale]),
            ),
        ];
        if !discrete_current_a.is_finite()
            || joules.iter().any(|v| !v.is_finite() || *v < 0.0)
            || accounts.iter().any(|&(residual, scale)| {
                let tolerance = self.controls.absolute_energy_tolerance_j
                    + self.controls.relative_balance_tolerance * scale;
                !residual.is_finite()
                    || !scale.is_finite()
                    || !tolerance.is_finite()
                    || residual.abs() > tolerance
            })
        {
            return Err(invalid(
                "nonlinear converter magnetic, bus, heat and work accounts do not close",
            )
            .into());
        }
        Ok(NonlinearConverterStep {
            native,
            magnetic_energy_change_j,
            bus_energy_change_j,
            terminal_work_j,
            conduction: ConductionHeat { joules },
            discrete_current_a,
        })
    }
}
