//! Bounded gate programs over one shared finite battery/motor native workspace.
use super::*;

impl PreparedCommutatedBatteryMotorCycle {
    /// Retain original source/lead, storage, motor and inertias. The explicit
    /// devices replace the supplied model's averaged modulation AND four bridge
    /// channel resistances; neither averaged channel heat nor ideal source work
    /// survives that replacement. Physical coordinates are never reinitialized.
    ///
    /// All constitutive laws are frozen over this cycle. At most 64 authored
    /// intervals and 1..=256 substeps per interval are admitted; each remainder
    /// must retain the original represented gate edge. Device support is checked
    /// at every substep's initial/discrete/final efforts, not just a cycle mean.
    pub fn new_controlled(
        model: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        if !(1..=MAX_SWITCH_INTERVALS).contains(&controls.max_intervals)
            || intervals.is_empty()
            || intervals.len() > controls.max_intervals
            || !(1..=MAX_POWER_SUBSTEPS_PER_INTERVAL).contains(&substeps_per_interval)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
            || !electrical.absolute_energy_tolerance_j.is_finite()
            || electrical.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&electrical.relative_energy_tolerance)
        {
            return Err(invalid(
                "commutated battery drive requires bounded intervals/refinement and finite energy controls",
            ).into());
        }
        electrical.solve.validate()?;
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .ok_or_else(|| invalid("commutated battery integration count overflow"))?;
        integration_steps
            .checked_mul(electrical.solve.max_iterations)
            .and_then(|n| n.checked_mul(usize::from(electrical.backtracking_halvings) + 1))
            .ok_or_else(|| invalid("commutated battery nonlinear work bound overflow"))?;
        let old = model.electrical().parameters();
        let parameters = CommutationParameters {
            inductance_h: old.inductance_h,
            bus_capacitance_f: old.bus_capacitance_f,
            winding_resistance_ohm: old.winding_resistance_ohm,
            devices,
        };
        parameters.validate()?;
        // This fixed native interconnection contains the original battery/lead,
        // converter winding, reciprocal motor and bearings. Bridge-only J/R are
        // zero here because the complete physical legs enter the supplied port.
        let battery = BatterySource::new(
            model.electrical().battery_curve().clone(),
            model.electrical().battery_parameters(),
        )?;
        let native = parameters
            .winding_model()?
            .with_battery_source(battery, model.electrical().lead_resistance_ohm())?
            .with_dc_motor(DcMotor::new(model.motor_parameters())?)?;
        let mut workspace = StepWorkspace::new(native.system())?;
        workspace.set_solve_budget(electrical.solve)?;
        workspace.set_backtracking_limit(electrical.backtracking_halvings)?;
        poll(&mut cancelled)?;
        let mut period = Dd::ZERO;
        let mut caches: Vec<CachedInterval> = Vec::new();
        let mut interval_steps = Vec::with_capacity(intervals.len());
        for interval in intervals {
            poll(&mut cancelled)?;
            if !positive(interval.duration_s) {
                return Err(
                    invalid("positive finite commutated battery gate duration required").into(),
                );
            }
            period = period + Dd::from_f64(interval.duration_s);
            let regular = interval.duration_s / substeps_per_interval as f64;
            let last = number(
                Dd::from_f64(interval.duration_s)
                    - Dd::from_f64(regular) * Dd::from_f64((substeps_per_interval - 1) as f64),
            )?;
            if !positive(regular) || !positive(last) {
                return Err(invalid("commutated battery subdivision is unrepresentable").into());
            }
            let mut represented = Dd::ZERO;
            for j in 0..substeps_per_interval {
                poll(&mut cancelled)?;
                represented = represented
                    + Dd::from_f64(if j + 1 == substeps_per_interval {
                        last
                    } else {
                        regular
                    });
            }
            if number(represented)?.to_bits() != interval.duration_s.to_bits() {
                return Err(invalid(
                    "commutated battery subdivisions must preserve each gate edge",
                )
                .into());
            }
            let mut indices = [0; 2];
            for (slot, duration_s) in [regular, last].into_iter().enumerate() {
                poll(&mut cancelled)?;
                let found = caches.iter().position(|c| {
                    c.gates == interval.gates && c.duration_s.to_bits() == duration_s.to_bits()
                });
                indices[slot] = if let Some(index) = found {
                    index
                } else {
                    let on = interval.gates.on;
                    let legs = [
                        HalfBridgeConduction::new(devices[0], devices[1], [on[0], on[1]])?,
                        HalfBridgeConduction::new(devices[2], devices[3], [on[2], on[3]])?,
                    ];
                    let index = caches.len();
                    caches.push(CachedInterval {
                        gates: interval.gates,
                        duration_s,
                        legs,
                    });
                    index
                };
            }
            interval_steps.push(IntervalSteps {
                regular: indices[0],
                last: indices[1],
            });
        }
        let period_s = number(period)?;
        if !positive(period_s) {
            return Err(invalid("commutated battery cycle duration is unrepresentable").into());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            model: native,
            parameters,
            workspace,
            intervals: intervals.to_vec(),
            interval_steps,
            caches,
            substeps_per_interval,
            integration_steps,
            period_s,
            electrical,
            controls,
        })
    }

    pub const fn parameters(&self) -> CommutationParameters {
        self.parameters
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.model.electrical().battery_parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.model.electrical().battery_curve()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.model.electrical().lead_resistance_ohm()
    }
    pub fn motor_parameters(&self) -> DcMotorParameters {
        self.model.motor_parameters()
    }
    pub fn intervals(&self) -> &[CommutationInterval] {
        &self.intervals
    }
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }
    pub const fn integration_steps(&self) -> usize {
        self.integration_steps
    }
    pub const fn substeps_per_interval(&self) -> usize {
        self.substeps_per_interval
    }
    /// Static leg/duration caches; all share ONE original nine-state workspace.
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    pub fn step_controlled(
        &mut self,
        initial: &[f64; 9],
        command: BatteryMotorSwitchCommand,
        next: &mut [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedBatteryMotorCycleStep, PreparedStepError> {
        self.advance(initial, |_| command, next, cancelled)
    }
    /// Commands belong to ORIGINAL intervals even when cached legs are shared.
    /// All nine states publish together after all native and whole-cycle gates.
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; 9],
        commands: &[BatteryMotorSwitchCommand],
        next: &mut [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedBatteryMotorCycleStep, PreparedStepError> {
        if commands.len() != self.intervals.len() {
            return Err(invalid(
                "one torque/load command is required per original commutation interval",
            )
            .into());
        }
        self.advance(initial, |k| commands[k], next, cancelled)
    }

    fn advance(
        &mut self,
        initial: &[f64; 9],
        command: impl Fn(usize) -> BatteryMotorSwitchCommand,
        next: &mut [f64; 9],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedBatteryMotorCycleStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let (mut changes, mut absolute_changes) = ([Dd::ZERO; 7], [Dd::ZERO; 7]);
        let (mut work, mut absolute_work) = ([Dd::ZERO; 8], [Dd::ZERO; 8]);
        let (mut means, mut native) = ([Dd::ZERO; 7], [Dd::ZERO; 3]);
        let mut heat = [Dd::ZERO; 14];
        let (mut curtailed, mut impulse) = (Dd::ZERO, Dd::ZERO);
        let (mut newton_iters, mut solver_residual) = (0usize, 0.0_f64);
        let mut last_means = [0.0; 7];
        let mut last_points = None;
        for k in 0..self.intervals.len() {
            poll(&mut cancelled)?;
            let command = command(k);
            let indices = self.interval_steps[k];
            for j in 0..self.substeps_per_interval {
                poll(&mut cancelled)?;
                let index = if j + 1 == self.substeps_per_interval {
                    indices.last
                } else {
                    indices.regular
                };
                let cache = &self.caches[index];
                let (endpoint, step) = Self::step_native(
                    &self.model,
                    self.parameters,
                    &mut self.workspace,
                    cache,
                    self.electrical,
                    &state,
                    command,
                    &mut cancelled,
                )?;
                let dt = Dd::from_f64(cache.duration_s);
                for i in 0..7 {
                    changes[i] = changes[i] + Dd::from_f64(step.storage_change_j[i]);
                    absolute_changes[i] =
                        absolute_changes[i] + Dd::from_f64(step.storage_change_j[i].abs());
                }
                for (sum, value) in heat.iter_mut().zip(step.heat_j) {
                    *sum = *sum + Dd::from_f64(value);
                }
                for (i, value) in step.work().into_iter().enumerate() {
                    work[i] = work[i] + Dd::from_f64(value);
                    absolute_work[i] = absolute_work[i] + Dd::from_f64(value.abs());
                }
                last_means = step.means();
                for i in 0..7 {
                    means[i] = means[i] + dt * Dd::from_f64(last_means[i]);
                }
                curtailed = curtailed + Dd::from_f64(step.curtailed_load_energy_j);
                impulse = impulse
                    + dt * (Dd::from_f64(command.shaft_torque_nm[0])
                        + Dd::from_f64(command.shaft_torque_nm[1]));
                for (sum, value) in native.iter_mut().zip([
                    step.native.delta_h,
                    step.native.dissipated,
                    step.native.supplied,
                ]) {
                    *sum = *sum + Dd::from_f64(value);
                }
                newton_iters = newton_iters
                    .checked_add(step.native.newton_iters)
                    .ok_or_else(|| invalid("commutated battery Newton count overflow"))?;
                solver_residual = solver_residual.max(step.native.solver_residual);
                last_points = Some((step.last_auxiliary_load, step.last_devices));
                state = endpoint;
            }
        }
        let p = self.parameters;
        let b = self.battery_parameters();
        let m = self.motor_parameters();
        let quadratic = |i, coefficient| quadratic_change(initial[i], state[i], coefficient);
        let storage_change_j = [
            self.battery_curve().energy_change_j(initial[0], state[0])?,
            quadratic(1, b.polarization_capacitance_f),
            quadratic(2, b.bus_capacitance_f),
            total([
                quadratic(3, p.inductance_h),
                quadratic(4, p.bus_capacitance_f[0]),
                quadratic(5, p.bus_capacitance_f[1]),
            ]),
            quadratic(6, m.inductance_h),
            quadratic(7, m.rotor_inertia_kg_m2),
            quadratic(8, m.stator_inertia_kg_m2),
        ];
        let delta_h = storage_energy_change(
            self.model.system().storage(),
            initial,
            &state,
            &mut [0.0; 9],
            &mut [0.0; 9],
        );
        let mean = if self.integration_steps == 1 {
            last_means
        } else {
            means.map(|v| (v / Dd::from_f64(self.period_s)).to_f64())
        };
        let work = work.map(Dd::to_f64);
        let absolute_work = absolute_work.map(Dd::to_f64);
        let absolute_changes = absolute_changes.map(Dd::to_f64);
        let (last_auxiliary_load, last_devices) =
            last_points.ok_or_else(|| invalid("empty commutated battery cycle"))?;
        let step = CommutatedBatteryMotorCycleStep {
            native: PreparedStepRecord {
                delta_h,
                dissipated: number(native[1])?,
                supplied: number(native[2])?,
                newton_iters,
                solver_residual,
                solver_certificate: None,
            },
            storage_change_j,
            battery_terminal_work_j: work[0],
            converter_terminal_work_j: [work[1], work[2]],
            motor_electrical_work_j: work[3],
            electromagnetic_work_j: work[4],
            shaft_work_j: [work[5], work[6]],
            load_work_j: work[7],
            curtailed_load_energy_j: number(curtailed)?,
            mean_bus_voltage_v: [mean[0], mean[1], mean[2]],
            mean_lead_current_a: mean[3],
            mean_motor_current_a: mean[4],
            mean_shaft_speed_rad_s: [mean[5], mean[6]],
            last_auxiliary_load,
            last_devices,
            angular_momentum_residual_nms: (Dd::from_f64(state[7]) - Dd::from_f64(initial[7])
                + Dd::from_f64(state[8])
                - Dd::from_f64(initial[8])
                - impulse)
                .to_f64(),
            integration_steps: self.integration_steps,
            heat_j: heat.map(Dd::to_f64),
        };
        step.admit(
            absolute_changes,
            absolute_work,
            self.controls.absolute_energy_tolerance_j,
            self.controls.relative_energy_tolerance,
        )?;
        for i in 0..7 {
            let residual = number(Dd::from_f64(storage_change_j[i]) - changes[i])?;
            let allowance = self.controls.absolute_energy_tolerance_j
                + self.controls.relative_energy_tolerance * absolute_changes[i];
            if !allowance.is_finite() || residual.abs() > allowance {
                return Err(invalid(
                    "commutated battery endpoint and accumulated constituent storage disagree",
                )
                .into());
            }
        }
        let scale = total(
            absolute_changes
                .into_iter()
                .chain(absolute_work)
                .chain(step.heat_j),
        );
        let allowance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residual = number(Dd::from_f64(delta_h) - native[0])?;
        if !scale.is_finite()
            || !allowance.is_finite()
            || residual.abs() > allowance
            || state.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "commutated battery endpoint and accumulated native energy disagree",
            )
            .into());
        }
        poll(&mut cancelled)?;
        *next = state;
        Ok(step)
    }
}
