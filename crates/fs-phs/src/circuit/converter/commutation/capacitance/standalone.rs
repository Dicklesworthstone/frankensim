//! Original five-/nine-state native commutation with finite-resistance terminals.
//! The gate lane appends four physical gate charges after grouped node charges;
//! its source work and four driver heats remain in the converter gate account.
use super::*;

#[derive(Clone, Debug)]
pub struct CapacitiveCommutatedStep {
    pub native: PreparedStepRecord,
    pub converter: CapacitiveConverterAccount,
    pub source_work_j: [f64; 2],
    pub terminal_resistor_heat_j: [f64; 2],
    pub load_work_j: [f64; 2],
    pub curtailed_load_energy_j: [f64; 2],
    pub mean_source_current_a: [f64; 2],
    pub last_loads: [PowerLoadPoint; 2],
    pub integration_steps: usize,
}
/// One reusable native workspace with explicitly continuing node charges.
/// Ideal supplies have finite impedance; they do not have finite stored energy.
pub struct PreparedCapacitiveCommutatedSwitchCycle<
    P = CapacitiveCommutationParameters,
    const N: usize = 5,
> {
    parameters: P,
    core: ConverterCore<N>,
    terminal_resistance_ohm: [f64; 2],
    plan: GatePlan,
    workspace: StepWorkspace,
}
struct Operating {
    nodes: NodePoint,
    loads: [PowerLoadPoint; 2],
    source: [f64; 2],
}
fn operating<const N: usize>(
    core: &ConverterCore<N>,
    gate: &GateStep,
    e: &[f64],
    sources: VoltageSources,
    loads: &[Option<RegulatedPowerLoad>; 2],
) -> Result<Operating, PhsError> {
    let nodes = core.point(gate, e)?;
    let loads = [load_point(loads[0], e[1])?, load_point(loads[1], e[2])?];
    let source = core::array::from_fn(|k| sources.current_a(k, e[k + 1]));
    if source.iter().any(|i| !i.is_finite()) {
        return Err(invalid("capacitive source current is unrepresentable"));
    }
    Ok(Operating {
        nodes,
        loads,
        source,
    })
}
/// The original charge curves select storage; gate and electrical step loops are shared.
pub type PreparedChargeCommutatedSwitchCycle =
    PreparedCapacitiveCommutatedSwitchCycle<charge::ChargeCommutationParameters>;

impl PreparedCapacitiveCommutatedSwitchCycle {
    pub fn new_controlled(
        parameters: CapacitiveCommutationParameters,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let core = ConverterCore::new(parameters)?;
        Self::prepare(
            parameters,
            core,
            terminal_resistance_ohm,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
impl PreparedCapacitiveCommutatedSwitchCycle<charge::ChargeCommutationParameters> {
    /// Retain four complete original V(q) curves. State coordinates must be
    /// initialized with these same curves; the constructor never converts,
    /// resets or silently clamps an existing constant-capacitor state.
    pub fn new_charge_controlled(
        parameters: charge::ChargeCommutationParameters,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let core = ConverterCore::new_charge(&parameters)?;
        Self::prepare(
            parameters,
            core,
            terminal_resistance_ohm,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
/// Four physical gate charges share the original converter trajectory. The
/// default channel is G(Vgs)*Vds; an immutable GateChannelResponse may instead
/// retain a full signed I(Vds,Vgs) source law with its analytic tangent.
pub type PreparedGateCommutatedSwitchCycle<C = GateConductanceCurve> =
    PreparedCapacitiveCommutatedSwitchCycle<gate::GateCommutationParameters<C>, 9>;
impl<C: GateChannelResponse + Clone + 'static>
    PreparedCapacitiveCommutatedSwitchCycle<gate::GateCommutationParameters<C>, 9>
{
    pub fn new_gate_controlled(
        parameters: gate::GateCommutationParameters<C>,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let core = ConverterCore::<9>::new_gate(&parameters)?;
        Self::prepare(
            parameters,
            core,
            terminal_resistance_ohm,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
impl<P: Clone, const N: usize> PreparedCapacitiveCommutatedSwitchCycle<P, N> {
    pub fn parameters(&self) -> P {
        self.parameters.clone()
    }
}
impl<P, const N: usize> PreparedCapacitiveCommutatedSwitchCycle<P, N> {
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        parameters: P,
        core: ConverterCore<N>,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        VoltageSources::new([0.0; 2], terminal_resistance_ohm)?;
        let plan = GatePlan::new(
            core.parameters.devices,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            core.gate_parameters().map(|p| &p.gate_drivers),
            &mut cancelled,
        )?;
        core.admit_work(&plan)?;
        let workspace = plan.workspace(&core.system)?;
        poll(&mut cancelled)?;
        Ok(Self {
            parameters,
            core,
            terminal_resistance_ohm,
            plan,
            workspace,
        })
    }
    pub fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.terminal_resistance_ohm
    }
    pub fn intervals(&self) -> &[CommutationInterval] {
        &self.plan.intervals
    }
    pub fn period_s(&self) -> f64 {
        self.plan.period_s
    }
    pub fn integration_steps(&self) -> usize {
        self.plan.integration_steps
    }
    pub fn substeps_per_interval(&self) -> usize {
        self.plan.substeps
    }
    pub fn cached_interval_count(&self) -> usize {
        self.plan.caches.len()
    }
    pub fn step_controlled(
        &mut self,
        initial: &[f64; N],
        source_voltage_v: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; N],
        mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveCommutatedStep, PreparedStepError> {
        let command = PowerSwitchCommand {
            source_voltage_v: *source_voltage_v,
            loads: *loads,
        };
        self.advance(initial, |_| command, next, mean_bus_voltage_v, cancelled)
    }
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; N],
        commands: &[PowerSwitchCommand],
        next: &mut [f64; N],
        mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveCommutatedStep, PreparedStepError> {
        if commands.len() != self.plan.intervals.len() {
            return Err(invalid(
                "one capacitive source/load command required per authored interval",
            )
            .into());
        }
        self.advance(
            initial,
            |k| commands[k],
            next,
            mean_bus_voltage_v,
            cancelled,
        )
    }
    fn advance(
        &mut self,
        initial: &[f64; N],
        command: impl Fn(usize) -> PowerSwitchCommand,
        next: &mut [f64; N],
        mean_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveCommutatedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let mut converter = ConverterAccumulator::new();
        let (
            mut source_work,
            mut resistor_heat,
            mut load_work,
            mut curtailed,
            mut current,
            mut transfers,
        ) = (
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
        );
        let mut native = [Dd::ZERO; 3];
        let mut gate_transfers = Dd::ZERO;
        let (mut iterations, mut residual) = (0usize, 0.0_f64);
        let mut last = None;
        for k in 0..self.plan.intervals.len() {
            poll(&mut cancelled)?;
            let command = command(k);
            let sources =
                VoltageSources::new(command.source_voltage_v, self.terminal_resistance_ohm)?;
            for j in 0..self.plan.substeps {
                poll(&mut cancelled)?;
                let gate = self.plan.step(k, j);
                operating(
                    &self.core,
                    gate,
                    &self.core.effort(&state)?,
                    sources,
                    &command.loads,
                )?;
                let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
                    let p = operating(&self.core, gate, e, sources, &command.loads).ok()?;
                    let heat = p.nodes.flow(out)?;
                    let mut supply = p.nodes.supplied_power();
                    for n in 0..2 {
                        let external =
                            Dd::from_f64(p.source[n]) - Dd::from_f64(p.loads[n].current_a);
                        out[1 + n] = number(p.nodes.legs[n].high_current - external).ok()?;
                        supply = supply + Dd::from_f64(e[1 + n]) * external;
                    }
                    Some(NonlinearPortPower {
                        supplied_power_w: number(supply).ok()?,
                        dissipated_power_w: heat,
                    })
                };
                let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
                    let Ok(p) = operating(&self.core, gate, e, sources, &command.loads) else {
                        return false;
                    };
                    if !p.nodes.tangent(de, out) {
                        return false;
                    }
                    for n in 0..2 {
                        out[n + 1] = total([
                            out[n + 1],
                            (sources.conductance_s[n] + p.loads[n].incremental_conductance_s)
                                * de[n + 1],
                        ]);
                    }
                    out.iter().all(|v| v.is_finite())
                };
                let (mut endpoint, mut outputs) = ([0.0; N], [0.0; 2]);
                let record = self.workspace.step_into_supplied_controlled(
                    &self.core.system,
                    &state,
                    &[0.0; 2],
                    gate.duration_s,
                    &mut endpoint,
                    &mut outputs,
                    &port,
                    None,
                    Some(&tangent),
                    &mut cancelled,
                )?;
                let e: &[f64; N] = self
                    .workspace
                    .last_step_discrete_effort()
                    .and_then(|e| e.try_into().ok())
                    .ok_or_else(|| {
                        invalid("capacitive step requires its original accepted effort")
                    })?;
                let p = operating(&self.core, gate, e, sources, &command.loads)?;
                operating(
                    &self.core,
                    gate,
                    &self.core.effort(&endpoint)?,
                    sources,
                    &command.loads,
                )?;
                let net = core::array::from_fn(|n| {
                    (Dd::from_f64(p.source[n]) - Dd::from_f64(p.loads[n].current_a)).to_f64()
                });
                let account = self.core.account(
                    &state,
                    &endpoint,
                    e,
                    p.nodes,
                    net,
                    gate.duration_s,
                    self.plan.electrical,
                )?;
                let dt = Dd::from_f64(gate.duration_s);
                let mut sw = [0.0; 2];
                let mut rh = [0.0; 2];
                let mut lw = [0.0; 2];
                for n in 0..2 {
                    sw[n] = product(gate.duration_s, command.source_voltage_v[n], p.source[n]);
                    rh[n] = (dt
                        * Dd::from_f64(self.terminal_resistance_ohm[n])
                        * Dd::from_f64(p.source[n])
                        * Dd::from_f64(p.source[n]))
                    .to_f64();
                    lw[n] = (dt * Dd::from_f64(p.loads[n].delivered_power_w)).to_f64();
                    let cu = (dt * Dd::from_f64(p.loads[n].curtailed_power_w)).to_f64();
                    if [rh[n], lw[n], cu]
                        .iter()
                        .any(|v| !v.is_finite() || *v < 0.0)
                    {
                        return Err(invalid(
                            "capacitive source/load heat and work are unrepresentable",
                        )
                        .into());
                    }
                    let scale =
                        total([sw[n].abs(), rh[n], lw[n], account.terminal_work_j[n].abs()]);
                    admit(
                        total([sw[n], -rh[n], -lw[n], -account.terminal_work_j[n]]),
                        scale,
                        self.plan.electrical.absolute_energy_tolerance_j,
                        self.plan.electrical.relative_energy_tolerance,
                        "capacitive source/resistor/load partition does not close",
                    )?;
                    source_work[n] = source_work[n] + Dd::from_f64(sw[n]);
                    resistor_heat[n] = resistor_heat[n] + Dd::from_f64(rh[n]);
                    load_work[n] = load_work[n] + Dd::from_f64(lw[n]);
                    curtailed[n] = curtailed[n] + Dd::from_f64(cu);
                    current[n] = current[n] + dt * Dd::from_f64(p.source[n]);
                    transfers[n] = transfers[n] + Dd::from_f64(scale);
                }
                gate_transfers = gate_transfers
                    + account
                        .gate_source_work_j()
                        .into_iter()
                        .fold(Dd::ZERO, |s, w| s + Dd::from_f64(w.abs()));
                let heat = total(account.all_heat_joules());
                let scale = total(sw.map(f64::abs))
                    .max(heat)
                    .max(record.delta_h.abs())
                    .max(total(account.gate_source_work_j().map(f64::abs)));
                for error in [
                    total([record.delta_h, -account.storage_change_j]),
                    total([record.dissipated, -heat]),
                    total(
                        account
                            .terminal_work_j
                            .into_iter()
                            .chain(account.gate_source_work_j())
                            .chain([-record.supplied]),
                    ),
                ] {
                    admit(
                        error,
                        scale,
                        self.plan.electrical.absolute_energy_tolerance_j,
                        self.plan.electrical.relative_energy_tolerance,
                        "native capacitive receipt and physical converter account disagree",
                    )?;
                }
                converter.push(&account, gate.duration_s);
                for (s, v) in
                    native
                        .iter_mut()
                        .zip([record.delta_h, record.dissipated, record.supplied])
                {
                    *s = *s + Dd::from_f64(v);
                }
                iterations = iterations
                    .checked_add(record.newton_iters)
                    .ok_or_else(|| invalid("capacitive Newton count overflow"))?;
                if !record.solver_residual.is_finite() {
                    return Err(invalid("nonfinite capacitive solver residual").into());
                }
                residual = residual.max(record.solver_residual);
                last = Some((p.loads, p.source));
                state = endpoint;
            }
        }
        let converter = converter.finish(&self.core, initial, &state, &self.plan)?;
        let summed_native_delta = number(native[0])?;
        let native = PreparedStepRecord {
            delta_h: converter.storage_change_j,
            dissipated: number(native[1])?,
            supplied: number(native[2])?,
            newton_iters: iterations,
            solver_residual: residual,
            solver_certificate: None,
        };
        let (last_loads, last_current) = last.ok_or_else(|| invalid("empty capacitive cycle"))?;
        let result = CapacitiveCommutatedStep {
            native,
            converter,
            source_work_j: source_work.map(Dd::to_f64),
            terminal_resistor_heat_j: resistor_heat.map(Dd::to_f64),
            load_work_j: load_work.map(Dd::to_f64),
            curtailed_load_energy_j: curtailed.map(Dd::to_f64),
            mean_source_current_a: if self.plan.integration_steps == 1 {
                last_current
            } else {
                current.map(|q| (q / Dd::from_f64(self.plan.period_s)).to_f64())
            },
            last_loads,
            integration_steps: self.plan.integration_steps,
        };
        let scale = number(transfers[0] + transfers[1] + gate_transfers)?
            .max(total(result.converter.all_heat_joules()));
        let c = self.plan.controls;
        for error in [
            total([result.native.delta_h, -summed_native_delta]),
            total([
                result.native.dissipated,
                -total(result.converter.all_heat_joules()),
            ]),
            total(
                result
                    .converter
                    .terminal_work_j
                    .into_iter()
                    .chain(result.converter.gate_source_work_j())
                    .chain([-result.native.supplied]),
            ),
            total(
                [result.native.delta_h, result.native.dissipated]
                    .into_iter()
                    .chain(result.terminal_resistor_heat_j)
                    .chain(result.load_work_j)
                    .chain(result.source_work_j.map(|v| -v))
                    .chain(result.converter.gate_source_work_j().map(|v| -v)),
            ),
        ] {
            admit(
                error,
                scale,
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "complete capacitive cycle native/source/storage balance does not close",
            )?;
        }
        for n in 0..2 {
            admit(
                total([
                    result.source_work_j[n],
                    -result.terminal_resistor_heat_j[n],
                    -result.load_work_j[n],
                    -result.converter.terminal_work_j[n],
                ]),
                number(transfers[n])?,
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "complete capacitive terminal partition does not close",
            )?;
        }
        if result
            .mean_source_current_a
            .iter()
            .chain(&result.curtailed_load_energy_j)
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("capacitive cycle mean or curtailed work is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        *next = state;
        *mean_bus_voltage_v = result.converter.mean_bus_voltage_v;
        Ok(result)
    }
}
