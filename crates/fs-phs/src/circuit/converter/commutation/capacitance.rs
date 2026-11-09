//! Constant output-capacitance charge storage in the original commutating bridge.
//!
//! State is [lambda, q_bus_in, q_bus_out, q_node_in, q_node_out].
//! A bus coordinate is TOTAL electrode charge Cb*Vbus+Ch*(Vbus-Vnode),
//! not the older bare-bus charge Cb*Vbus. Node charge is
//! -Ch*Vbus+(Ch+Cl)*Vnode. Initial node voltages must be authored explicitly;
//! gate edges never reset charge or voltage, and no Coss energy is debited as
//! an additional switching-loss map.
//!
//! Each grounded reciprocal block is admitted and inverted by the existing
//! CapacitanceBlock owner. The same original quadratic storage and native
//! solver serve both the standalone five-state and finite-source eleven-state
//! lanes. Constant positive device capacitances are required: nonlinear C(V),
//! temperature-dependent storage, recovery charge and gate-drive energy are
//! not inferred. Small RC time scales require explicit temporal refinement.
//!
//! The `charge` module admits complete original nonlinear device Q(V) laws
//! and their coupled electrode storage without reducing them to constant C.

use super::{
    CommutationInterval, CommutationParameters, DeviceConductionPoint, DeviceGateState,
    HalfBridgeConduction, HalfBridgeNodePoint, SwitchDevice, invalid, number, poll, positive,
    total,
};
use crate::circuit::capacitance::CapacitanceBlock;
use crate::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadPoint, RegulatedPowerLoad, VoltageSources, load_point,
};
use crate::circuit::converter::switching::power_load::{
    MAX_POWER_SUBSTEPS_PER_INTERVAL, PowerSwitchCommand,
};
use crate::circuit::converter::switching::{MAX_SWITCH_INTERVALS, SwitchCycleControls};
use crate::{
    NonlinearPortPower, PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord,
    QuadraticStorage, StepWorkspace, Storage, storage_energy_change,
};
use fs_math::dd::Dd;

/// Original finite battery and motor shafts with continuing device capacitor charge.
pub mod battery_motor;
/// Original voltage-dependent device charge and exact electrode storage.
pub mod charge;
mod plan;
mod standalone;
use plan::{GatePlan, GateStep};
pub use standalone::{CapacitiveCommutatedStep, PreparedCapacitiveCommutatedSwitchCycle};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapacitiveCommutationParameters {
    pub converter: CommutationParameters,
    /// Constant positive C in input-high/input-low/output-high/output-low order.
    /// These are physical drain/source storage, not differential C(V) samples.
    pub output_capacitance_f: [f64; 4],
}
impl CapacitiveCommutationParameters {
    pub fn validate(&self) -> Result<(), PhsError> {
        ConverterCore::new(*self).map(|_| ())
    }
    /// Construct total electrode charges from explicitly authored physical
    /// voltages. Old three-state/drive bus charges must not be copied unchanged
    /// into these coordinates: the high-device capacitor shares the bus node.
    pub fn state(
        &self,
        current_a: f64,
        bus_voltage_v: [f64; 2],
        node_voltage_v: [f64; 2],
    ) -> Result<[f64; 5], PhsError> {
        ConverterCore::new(*self)?.state(current_a, bus_voltage_v, node_voltage_v)
    }
    /// Cold source-binding query using the SAME original coupled storage law.
    /// Gate-specific device support is checked when the authored program steps.
    pub fn stored_energy_j(&self, state: &[f64; 5]) -> Result<f64, PhsError> {
        ConverterCore::new(*self)?.energy(state)
    }
}

/// Shared converter account on either the standalone or finite-source trajectory.
/// Capacitor energy changes are reversible storage, never extra heat channels.
#[derive(Clone, Debug)]
pub struct CapacitiveConverterAccount {
    pub storage_change_j: f64,
    pub inductor_energy_change_j: f64,
    pub bus_capacitor_energy_change_j: [f64; 2],
    pub device_capacitor_energy_change_j: [f64; 4],
    pub terminal_work_j: [f64; 2],
    pub winding_heat_j: f64,
    pub channel_heat_j: [f64; 4],
    pub diode_heat_j: [f64; 4],
    pub mean_current_a: f64,
    pub mean_bus_voltage_v: [f64; 2],
    pub mean_node_voltage_v: [f64; 2],
    pub last_devices: [DeviceConductionPoint; 4],
}
impl CapacitiveConverterAccount {
    /// Winding, four channels, four diodes. Coss is absent because it is storage.
    pub fn heat_joules(&self) -> [f64; 9] {
        let c = self.channel_heat_j;
        let d = self.diode_heat_j;
        [
            self.winding_heat_j,
            c[0],
            c[1],
            c[2],
            c[3],
            d[0],
            d[1],
            d[2],
            d[3],
        ]
    }
    pub fn balance_residual_j(&self) -> f64 {
        total(self.heat_joules().into_iter().chain([
            self.storage_change_j,
            -self.terminal_work_j[0],
            -self.terminal_work_j[1],
        ]))
    }
    fn storage_parts(&self) -> [f64; 7] {
        let b = self.bus_capacitor_energy_change_j;
        let d = self.device_capacitor_energy_change_j;
        [
            self.inductor_energy_change_j,
            b[0],
            b[1],
            d[0],
            d[1],
            d[2],
            d[3],
        ]
    }
    fn means(&self) -> [f64; 5] {
        [
            self.mean_current_a,
            self.mean_bus_voltage_v[0],
            self.mean_bus_voltage_v[1],
            self.mean_node_voltage_v[0],
            self.mean_node_voltage_v[1],
        ]
    }
    fn admit(&self, scale: f64, absolute: f64, relative: f64) -> Result<(), PhsError> {
        let parts = self.storage_parts();
        let heat = self.heat_joules();
        if parts
            .iter()
            .chain(&self.terminal_work_j)
            .chain(&self.means())
            .any(|v| !v.is_finite())
            || heat.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "capacitive converter work, energy or heat is unrepresentable",
            ));
        }
        let scale = scale.max(total(parts.map(f64::abs))).max(total(heat));
        admit(
            self.balance_residual_j(),
            scale,
            absolute,
            relative,
            "capacitive converter storage, terminal work and physical heat do not close",
        )?;
        admit(
            total(parts.into_iter().chain([-self.storage_change_j])),
            scale,
            absolute,
            relative,
            "original bus/device capacitor and winding storage partition does not close",
        )
    }
}
fn admit(
    residual: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
    what: &'static str,
) -> Result<(), PhsError> {
    let allowance = absolute + relative * scale;
    if !residual.is_finite()
        || !scale.is_finite()
        || !allowance.is_finite()
        || residual.abs() > allowance
    {
        Err(invalid(what))
    } else {
        Ok(())
    }
}
fn product(a: f64, b: f64, c: f64) -> f64 {
    (Dd::from_f64(a) * Dd::from_f64(b) * Dd::from_f64(c)).to_f64()
}
fn quadratic(before: f64, after: f64, coefficient: f64) -> f64 {
    let a = Dd::from_f64(before);
    let b = Dd::from_f64(after);
    ((b - a) * (Dd::from_f64(0.5) * a + Dd::from_f64(0.5) * b) / Dd::from_f64(coefficient)).to_f64()
}
fn capacitor_change(before: Dd, after: Dd, capacitance: f64) -> f64 {
    (Dd::from_f64(capacitance)
        * (after - before)
        * (Dd::from_f64(0.5) * before + Dd::from_f64(0.5) * after))
        .to_f64()
}

struct ConverterCore {
    parameters: CapacitiveCommutationParameters,
    blocks: [CapacitanceBlock; 2],
    system: PortHamiltonian,
}
impl ConverterCore {
    fn new(parameters: CapacitiveCommutationParameters) -> Result<Self, PhsError> {
        let p = parameters.converter;
        p.validate()?;
        if parameters
            .output_capacitance_f
            .iter()
            .any(|c| !positive(*c) || !positive(1.0 / c))
        {
            return Err(invalid(
                "four authored positive representable constant output capacitances required",
            ));
        }
        let block = |leg: usize| {
            let high = parameters.output_capacitance_f[2 * leg];
            let low = parameters.output_capacitance_f[2 * leg + 1];
            let bus = p.bus_capacitance_f[leg];
            CapacitanceBlock::try_new(
                vec![1 + leg, 3 + leg],
                vec![bus + high, -high, -high, high + low],
            )
            .map_err(|_| {
                invalid("grounded output/bus capacitance block is not representable SPD storage")
            })
        };
        let blocks = [block(0)?, block(1)?];
        let mut q = vec![0.0; 25];
        q[0] = 1.0 / p.inductance_h;
        for (leg, block) in blocks.iter().enumerate() {
            let indices = [1 + leg, 3 + leg];
            let inverse = block.inverse_matrix_per_f();
            for i in 0..2 {
                for j in 0..2 {
                    q[5 * indices[i] + indices[j]] = inverse[2 * i + j];
                }
            }
        }
        // Direct sum of the two already admitted reciprocal SPD inverses and
        // positive inductor compliance. No additional inverse or fitted C law.
        let storage = QuadraticStorage::admitted(q, 5);
        let mut j = vec![0.0; 25];
        j[3] = 1.0;
        j[4] = -1.0;
        j[15] = -1.0;
        j[20] = 1.0;
        let mut r = vec![0.0; 25];
        r[0] = p.winding_resistance_ohm;
        let mut g = vec![0.0; 10];
        g[2] = 1.0;
        g[5] = 1.0;
        let system = PortHamiltonian::new(5, 2, j, r, g, Box::new(storage))?;
        Ok(Self {
            parameters,
            blocks,
            system,
        })
    }
    fn state(&self, current: f64, bus: [f64; 2], node: [f64; 2]) -> Result<[f64; 5], PhsError> {
        if !current.is_finite() || bus.iter().chain(&node).any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite authored current and bus/node voltages required",
            ));
        }
        let charges = |leg: usize| {
            self.blocks[leg]
                .charges_c(&[bus[leg], node[leg]])
                .map_err(|_| invalid("original bus/node electrode charges are unrepresentable"))
        };
        let a = charges(0)?;
        let b = charges(1)?;
        let state = [
            self.parameters.converter.inductance_h * current,
            a[0],
            b[0],
            a[1],
            b[1],
        ];
        self.energy(&state)?;
        Ok(state)
    }
    fn energy(&self, state: &[f64; 5]) -> Result<f64, PhsError> {
        if state.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite original capacitive converter state required",
            ));
        }
        let h = self.system.hamiltonian(state);
        if !h.is_finite() || h < 0.0 || (h == 0.0 && state.iter().any(|v| *v != 0.0)) {
            return Err(invalid(
                "nonzero original capacitive converter energy is unrepresentable",
            ));
        }
        Ok(h)
    }
    fn effort(&self, state: &[f64; 5]) -> Result<[f64; 5], PhsError> {
        self.energy(state)?;
        let mut e = [0.0; 5];
        self.system.storage().gradient(state, &mut e);
        if e.iter().any(|v| !v.is_finite()) {
            return Err(invalid("capacitive converter effort is unrepresentable"));
        }
        Ok(e)
    }
    fn storage_parts(
        &self,
        initial: &[f64; 5],
        next: &[f64; 5],
    ) -> Result<([f64; 7], f64), PhsError> {
        let a = self.effort(initial)?;
        let b = self.effort(next)?;
        let p = self.parameters;
        let mut parts = [0.0; 7];
        parts[0] = quadratic(initial[0], next[0], p.converter.inductance_h);
        for k in 0..2 {
            let av = Dd::from_f64(a[1 + k]);
            let bv = Dd::from_f64(b[1 + k]);
            let an = Dd::from_f64(a[3 + k]);
            let bn = Dd::from_f64(b[3 + k]);
            parts[1 + k] = capacitor_change(av, bv, p.converter.bus_capacitance_f[k]);
            parts[3 + 2 * k] = capacitor_change(av - an, bv - bn, p.output_capacitance_f[2 * k]);
            parts[4 + 2 * k] = capacitor_change(an, bn, p.output_capacitance_f[2 * k + 1]);
        }
        let delta = storage_energy_change(
            self.system.storage(),
            initial,
            next,
            &mut [0.0; 5],
            &mut [0.0; 5],
        );
        Ok((parts, delta))
    }
    fn account(
        &self,
        initial: &[f64; 5],
        next: &[f64; 5],
        e: &[f64; 5],
        point: NodePoint,
        current: [f64; 2],
        dt: f64,
        controls: PowerLoadControls,
    ) -> Result<CapacitiveConverterAccount, PhsError> {
        let (s, delta) = self.storage_parts(initial, next)?;
        let d = point.devices();
        let account = CapacitiveConverterAccount {
            storage_change_j: delta,
            inductor_energy_change_j: s[0],
            bus_capacitor_energy_change_j: [s[1], s[2]],
            device_capacitor_energy_change_j: [s[3], s[4], s[5], s[6]],
            terminal_work_j: [product(dt, e[1], current[0]), product(dt, e[2], current[1])],
            winding_heat_j: (Dd::from_f64(dt)
                * Dd::from_f64(self.parameters.converter.winding_resistance_ohm)
                * Dd::from_f64(e[0])
                * Dd::from_f64(e[0]))
            .to_f64(),
            channel_heat_j: d
                .map(|d| (Dd::from_f64(dt) * Dd::from_f64(d.channel_power_w)).to_f64()),
            diode_heat_j: d.map(|d| (Dd::from_f64(dt) * Dd::from_f64(d.diode_power_w)).to_f64()),
            mean_current_a: e[0],
            mean_bus_voltage_v: [e[1], e[2]],
            mean_node_voltage_v: [e[3], e[4]],
            last_devices: d,
        };
        account.admit(
            total(account.terminal_work_j.map(f64::abs)),
            controls.absolute_energy_tolerance_j,
            controls.relative_energy_tolerance,
        )?;
        Ok(account)
    }
}

#[derive(Clone, Copy)]
struct NodePoint {
    legs: [HalfBridgeNodePoint; 2],
}
impl NodePoint {
    fn devices(self) -> [DeviceConductionPoint; 4] {
        [
            self.legs[0].high,
            self.legs[0].low,
            self.legs[1].high,
            self.legs[1].low,
        ]
    }
    fn flow(self, out: &mut [f64]) -> Option<f64> {
        if out.len() != 5 {
            return None;
        }
        out.fill(0.0);
        for k in 0..2 {
            out[1 + k] = number(self.legs[k].high_current).ok()?;
            out[3 + k] = number(self.legs[k].low_current - self.legs[k].high_current).ok()?;
        }
        number(self.legs[0].loss + self.legs[1].loss).ok()
    }
    fn tangent(self, de: &[f64], out: &mut [f64]) -> bool {
        if de.len() != 5 || out.len() != 5 {
            return false;
        }
        out.fill(0.0);
        for k in 0..2 {
            let high = Dd::from_f64(self.legs[k].d_high_current_d_drain_source_voltage_s)
                * (Dd::from_f64(de[1 + k]) - Dd::from_f64(de[3 + k]));
            let low = Dd::from_f64(self.legs[k].d_low_current_d_drain_source_voltage_s)
                * Dd::from_f64(de[3 + k]);
            out[1 + k] = high.to_f64();
            out[3 + k] = (low - high).to_f64();
        }
        out.iter().all(|v| v.is_finite())
    }
}
impl GateStep {
    fn point(&self, e: &[f64]) -> Result<NodePoint, PhsError> {
        if e.len() != 5 || e.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite original capacitive converter efforts required",
            ));
        }
        Ok(NodePoint {
            legs: [
                self.legs[0].evaluate_at_node(e[1], e[3])?,
                self.legs[1].evaluate_at_node(e[2], e[4])?,
            ],
        })
    }
}

struct ConverterAccumulator {
    heat: [Dd; 9],
    work: [Dd; 2],
    absolute_work: Dd,
    parts: [Dd; 7],
    absolute_parts: [Dd; 7],
    delta: Dd,
    means: [Dd; 5],
    last: Option<([f64; 5], [DeviceConductionPoint; 4])>,
}
impl ConverterAccumulator {
    fn new() -> Self {
        Self {
            heat: [Dd::ZERO; 9],
            work: [Dd::ZERO; 2],
            absolute_work: Dd::ZERO,
            parts: [Dd::ZERO; 7],
            absolute_parts: [Dd::ZERO; 7],
            delta: Dd::ZERO,
            means: [Dd::ZERO; 5],
            last: None,
        }
    }
    fn push(&mut self, a: &CapacitiveConverterAccount, dt: f64) {
        for (s, v) in self.heat.iter_mut().zip(a.heat_joules()) {
            *s = *s + Dd::from_f64(v);
        }
        for k in 0..2 {
            self.work[k] = self.work[k] + Dd::from_f64(a.terminal_work_j[k]);
            self.absolute_work = self.absolute_work + Dd::from_f64(a.terminal_work_j[k].abs());
        }
        for (k, v) in a.storage_parts().into_iter().enumerate() {
            self.parts[k] = self.parts[k] + Dd::from_f64(v);
            self.absolute_parts[k] = self.absolute_parts[k] + Dd::from_f64(v.abs());
        }
        self.delta = self.delta + Dd::from_f64(a.storage_change_j);
        let means = a.means();
        for k in 0..5 {
            self.means[k] = self.means[k] + Dd::from_f64(dt) * Dd::from_f64(means[k]);
        }
        self.last = Some((means, a.last_devices));
    }
    fn finish(
        self,
        core: &ConverterCore,
        initial: &[f64; 5],
        next: &[f64; 5],
        plan: &GatePlan,
    ) -> Result<CapacitiveConverterAccount, PhsError> {
        let (s, delta) = core.storage_parts(initial, next)?;
        let (last, last_devices) = self
            .last
            .ok_or_else(|| invalid("empty capacitor-account cycle"))?;
        let means = if plan.integration_steps == 1 {
            last
        } else {
            self.means
                .map(|v| (v / Dd::from_f64(plan.period_s)).to_f64())
        };
        let heat = self.heat.map(Dd::to_f64);
        let a = CapacitiveConverterAccount {
            storage_change_j: delta,
            inductor_energy_change_j: s[0],
            bus_capacitor_energy_change_j: [s[1], s[2]],
            device_capacitor_energy_change_j: [s[3], s[4], s[5], s[6]],
            terminal_work_j: self.work.map(Dd::to_f64),
            winding_heat_j: heat[0],
            channel_heat_j: [heat[1], heat[2], heat[3], heat[4]],
            diode_heat_j: [heat[5], heat[6], heat[7], heat[8]],
            mean_current_a: means[0],
            mean_bus_voltage_v: [means[1], means[2]],
            mean_node_voltage_v: [means[3], means[4]],
            last_devices,
        };
        let scale = number(
            self.absolute_work + self.absolute_parts.into_iter().fold(Dd::ZERO, |a, b| a + b),
        )?;
        let c = plan.controls;
        a.admit(
            scale,
            c.absolute_energy_tolerance_j,
            c.relative_energy_tolerance,
        )?;
        for k in 0..7 {
            admit(
                number(Dd::from_f64(s[k]) - self.parts[k])?,
                number(self.absolute_parts[k])?,
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "capacitor endpoint and accumulated physical storage disagree",
            )?;
        }
        admit(
            number(Dd::from_f64(delta) - self.delta)?,
            scale,
            c.absolute_energy_tolerance_j,
            c.relative_energy_tolerance,
            "capacitive converter endpoint and accumulated native storage disagree",
        )?;
        Ok(a)
    }
}
