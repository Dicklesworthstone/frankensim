//! Physical output-capacitor charge in the original commutating bridge.
//!
//! The original five-state lanes retain [lambda, q_bus_in, q_bus_out,
//! q_node_in, q_node_out]. Bus charge is TOTAL electrode charge
//! Cb*Vbus+Qhigh(Vbus-Vnode); node charge is Qlow(Vnode)-Qhigh(Vbus-Vnode).
//! The gate-driven lane has nine states: flux, two bus charges, two GROUPED
//! node charges, and four gate charges. A grouped node includes its high-side
//! gate electrode, making each gate effort the actual source-referenced Vgs.
//! Initialize through the retained original parameters. Bare-bus, ungrouped
//! node, and differently authored capacitor-law charges are not interchangeable.
//!
//! Constant capacitors retain CapacitanceBlock; original Q(V) curves retain
//! their bounded source laws and AVF rule. The reciprocal gate wrapper adds
//! authored Cgs/Cgd through the same storage owner. Five/nine-state converters
//! and eleven/fifteen-state finite drives share native stepping and accounts.
//! Gate commands drive finite resistors; they do not instantaneously overwrite
//! gate charge or channel conductance. Actual gate supply work and driver heat
//! are explicit. The driven channel may retain a complete signed I(Vds,Vgs)
//! source response, including saturation and its original two-voltage tangent;
//! the existing gate-only conductance curves remain the default. Channel and
//! body-diode current are distinct physical owners. No Coss loss map, recovery charge, temperature-varying storage
//! or finite stored gate-rail energy is inferred. Small RC scales still require
//! explicitly authored temporal refinement.

use super::gate::{
    GateChannelLaw, GateChannelResponse, GateConductanceCurve, GateDriver, GateDriverBranch,
    GateDriverPoint,
};
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
/// Original bounded device charge curves and exact coupled electrode storage.
pub mod charge;
/// Reciprocal intrinsic drain, gate/source and Miller charge.
pub mod gate;
mod plan;
mod standalone;
pub use gate::GatePowerAccount;
use plan::{GatePlan, GateStep};
pub use standalone::{
    CapacitiveCommutatedStep, PreparedCapacitiveCommutatedSwitchCycle,
    PreparedChargeCommutatedSwitchCycle, PreparedGateCommutatedSwitchCycle,
};

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
    pub gate: Option<GatePowerAccount>,
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
    /// Complete original device/winding heat plus four physical gate drivers.
    pub fn all_heat_joules(&self) -> [f64; 13] {
        let mut out = [0.0; 13];
        out[..9].copy_from_slice(&self.heat_joules());
        if let Some(gate) = &self.gate {
            out[9..].copy_from_slice(&gate.resistor_heat_j);
        }
        out
    }
    pub fn gate_source_work_j(&self) -> [f64; 4] {
        self.gate
            .as_ref()
            .map(|g| g.source_work_j)
            .unwrap_or([0.0; 4])
    }
    pub fn balance_residual_j(&self) -> f64 {
        total(
            self.all_heat_joules()
                .into_iter()
                .chain([
                    self.storage_change_j,
                    -self.terminal_work_j[0],
                    -self.terminal_work_j[1],
                ])
                .chain(self.gate_source_work_j().map(|w| -w)),
        )
    }
    fn storage_parts(&self) -> [f64; 15] {
        let b = self.bus_capacitor_energy_change_j;
        let d = self.device_capacitor_energy_change_j;
        let mut out = [0.0; 15];
        out[..7].copy_from_slice(&[
            self.inductor_energy_change_j,
            b[0],
            b[1],
            d[0],
            d[1],
            d[2],
            d[3],
        ]);
        if let Some(gate) = &self.gate {
            out[7..11].copy_from_slice(&gate.gate_source_energy_change_j);
            out[11..].copy_from_slice(&gate.gate_drain_energy_change_j);
        }
        out
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
        let heat = self.all_heat_joules();
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
        if let Some(gate) = &self.gate {
            if gate
                .source_work_j
                .iter()
                .chain(&gate.mean_voltage_v)
                .chain(&gate.mean_current_a)
                .chain(&gate.last_driver_current_a)
                .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "gate source work or conjugate effort is unrepresentable",
                ));
            }
        }
        let scale = scale
            .max(total(parts.map(f64::abs)))
            .max(total(heat))
            .max(total(self.gate_source_work_j().map(f64::abs)));
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

enum ConverterStorage {
    Constant {
        output_capacitance_f: [f64; 4],
        blocks: [CapacitanceBlock; 2],
    },
    Charge(charge::ChargeCommutationStorage),
    Gate(gate::GateCommutationStorage),
}
struct ConverterCore<const N: usize = 5> {
    parameters: CommutationParameters,
    storage: ConverterStorage,
    system: PortHamiltonian,
}
impl ConverterCore<5> {
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
        // Existing admitted reciprocal SPD blocks remain the constant-law path.
        Self::with_storage(
            p,
            ConverterStorage::Constant {
                output_capacitance_f: parameters.output_capacitance_f,
                blocks,
            },
            Box::new(QuadraticStorage::admitted(q, 5)),
        )
    }
    fn new_charge(parameters: &charge::ChargeCommutationParameters) -> Result<Self, PhsError> {
        let storage = charge::ChargeCommutationStorage::new(parameters)?;
        Self::with_storage(
            parameters.converter,
            ConverterStorage::Charge(storage.clone()),
            Box::new(storage),
        )
    }
    fn state(&self, current: f64, bus: [f64; 2], node: [f64; 2]) -> Result<[f64; 5], PhsError> {
        if let ConverterStorage::Charge(storage) = &self.storage {
            return storage.state(current, bus, node);
        }
        let ConverterStorage::Constant { blocks, .. } = &self.storage else {
            unreachable!()
        };
        if !current.is_finite() || bus.iter().chain(&node).any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite authored current and bus/node voltages required",
            ));
        }
        let charges = |leg: usize| {
            blocks[leg]
                .charges_c(&[bus[leg], node[leg]])
                .map_err(|_| invalid("original bus/node electrode charges are unrepresentable"))
        };
        let a = charges(0)?;
        let b = charges(1)?;
        let state = [
            self.parameters.inductance_h * current,
            a[0],
            b[0],
            a[1],
            b[1],
        ];
        self.energy(&state)?;
        Ok(state)
    }
}
impl ConverterCore<9> {
    fn new_gate<C: GateChannelResponse + Clone + 'static>(
        parameters: &gate::GateCommutationParameters<C>,
    ) -> Result<Self, PhsError> {
        let storage = gate::GateCommutationStorage::new(parameters)?;
        Self::with_storage(
            parameters.converter,
            ConverterStorage::Gate(storage.clone()),
            Box::new(storage),
        )
    }
}
impl<const N: usize> ConverterCore<N> {
    fn with_storage(
        parameters: CommutationParameters,
        storage: ConverterStorage,
        law: Box<dyn Storage>,
    ) -> Result<Self, PhsError> {
        let expected = if matches!(&storage, ConverterStorage::Gate(_)) {
            9
        } else {
            5
        };
        if N != expected {
            return Err(invalid("original converter storage dimension mismatch"));
        }
        let mut j = vec![0.0; N * N];
        j[3] = 1.0;
        j[4] = -1.0;
        j[3 * N] = -1.0;
        j[4 * N] = 1.0;
        let mut r = vec![0.0; N * N];
        r[0] = parameters.winding_resistance_ohm;
        let mut g = vec![0.0; N * 2];
        g[2] = 1.0;
        g[5] = 1.0;
        let system = PortHamiltonian::new(N, 2, j, r, g, law)?;
        Ok(Self {
            parameters,
            storage,
            system,
        })
    }
    fn gate_parameters(&self) -> Option<&gate::GateCommutationParameters<GateChannelLaw>> {
        if let ConverterStorage::Gate(storage) = &self.storage {
            Some(storage.parameters())
        } else {
            None
        }
    }
    fn admit_work(&self, plan: &GatePlan) -> Result<(), PhsError> {
        let cells = match &self.storage {
            ConverterStorage::Charge(storage) => storage.maximum_cells_per_chord(),
            ConverterStorage::Gate(storage) => storage.maximum_cells_per_chord(),
            ConverterStorage::Constant { .. } => 1,
        };
        plan.integration_steps
            .checked_mul(plan.electrical.solve.max_iterations)
            .and_then(|n| n.checked_mul(usize::from(plan.electrical.backtracking_halvings) + 1))
            .and_then(|n| n.checked_mul(2))
            .and_then(|n| n.checked_mul(cells))
            .ok_or_else(|| invalid("original converter chord/cycle work bound overflow"))?;
        Ok(())
    }
    fn energy(&self, state: &[f64; N]) -> Result<f64, PhsError> {
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
    fn effort(&self, state: &[f64; N]) -> Result<[f64; N], PhsError> {
        self.energy(state)?;
        let mut e = [0.0; N];
        self.system.storage().gradient(state, &mut e);
        if e.iter().any(|v| !v.is_finite()) {
            return Err(invalid("capacitive converter effort is unrepresentable"));
        }
        Ok(e)
    }
    fn storage_parts(
        &self,
        initial: &[f64; N],
        next: &[f64; N],
    ) -> Result<([f64; 15], f64), PhsError> {
        let parts = match &self.storage {
            ConverterStorage::Charge(storage) => {
                let values = storage.parts(
                    initial
                        .as_slice()
                        .try_into()
                        .map_err(|_| invalid("charge dimension"))?,
                    next.as_slice()
                        .try_into()
                        .map_err(|_| invalid("charge dimension"))?,
                )?;
                let mut out = [0.0; 15];
                out[..7].copy_from_slice(&values);
                out
            }
            ConverterStorage::Gate(storage) => storage.parts(
                initial
                    .as_slice()
                    .try_into()
                    .map_err(|_| invalid("gate charge dimension"))?,
                next.as_slice()
                    .try_into()
                    .map_err(|_| invalid("gate charge dimension"))?,
            )?,
            ConverterStorage::Constant {
                output_capacitance_f,
                ..
            } => {
                let a = self.effort(initial)?;
                let b = self.effort(next)?;
                let p = self.parameters;
                let mut parts = [0.0; 15];
                parts[0] = quadratic(initial[0], next[0], p.inductance_h);
                for k in 0..2 {
                    let av = Dd::from_f64(a[1 + k]);
                    let bv = Dd::from_f64(b[1 + k]);
                    let an = Dd::from_f64(a[3 + k]);
                    let bn = Dd::from_f64(b[3 + k]);
                    parts[1 + k] = capacitor_change(av, bv, p.bus_capacitance_f[k]);
                    parts[3 + 2 * k] =
                        capacitor_change(av - an, bv - bn, output_capacitance_f[2 * k]);
                    parts[4 + 2 * k] = capacitor_change(an, bn, output_capacitance_f[2 * k + 1]);
                }
                parts
            }
        };
        let delta = storage_energy_change(
            self.system.storage(),
            initial,
            next,
            &mut [0.0; N],
            &mut [0.0; N],
        );
        Ok((parts, delta))
    }
    fn account(
        &self,
        initial: &[f64; N],
        next: &[f64; N],
        e: &[f64; N],
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
                * Dd::from_f64(self.parameters.winding_resistance_ohm)
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
            gate: point.gate.map(|g| GatePowerAccount {
                gate_source_energy_change_j: core::array::from_fn(|k| s[7 + k]),
                gate_drain_energy_change_j: core::array::from_fn(|k| s[11 + k]),
                source_work_j: g.driver.map(|p| (Dd::from_f64(dt) * p.supplied).to_f64()),
                resistor_heat_j: g.driver.map(|p| (Dd::from_f64(dt) * p.loss).to_f64()),
                mean_voltage_v: core::array::from_fn(|k| e[5 + k]),
                mean_current_a: g.driver.map(|p| p.current_a),
                last_driver_current_a: g.driver.map(|p| p.current_a),
            }),
        };
        account.admit(
            total(account.terminal_work_j.map(f64::abs)),
            controls.absolute_energy_tolerance_j,
            controls.relative_energy_tolerance,
        )?;
        Ok(account)
    }
    fn point(&self, gate: &GateStep, e: &[f64]) -> Result<NodePoint, PhsError> {
        let Some(p) = self.gate_parameters() else {
            return gate.point(e);
        };
        if e.len() != 9 || e.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite nine-state gate efforts required"));
        }
        let voltage = [
            number(Dd::from_f64(e[1]) - Dd::from_f64(e[3]))?,
            e[3],
            number(Dd::from_f64(e[2]) - Dd::from_f64(e[4]))?,
            e[4],
        ];
        let c = four(core::array::from_fn(|k| {
            p.channel_conductance[k].evaluate(voltage[k], e[5 + k])
        }))?;
        let branches = gate
            .drivers
            .as_ref()
            .ok_or_else(|| invalid("gate driver branches not prepared"))?;
        let driver = four(core::array::from_fn(|k| branches[k].evaluate(e[5 + k])))?;
        let legs = [
            gate.legs[0].evaluate_at_node_with_channels(e[1], e[3], [c[0], c[1]])?,
            gate.legs[1].evaluate_at_node_with_channels(e[2], e[4], [c[2], c[3]])?,
        ];
        Ok(NodePoint {
            legs,
            gate: Some(DrivenPoint {
                driver,
                channel_slope: c.map(|point| point.gate_derivative_s),
            }),
        })
    }
}

fn four<T>(values: [Result<T, PhsError>; 4]) -> Result<[T; 4], PhsError> {
    let [a, b, c, d] = values;
    Ok([a?, b?, c?, d?])
}
#[derive(Clone, Copy)]
struct DrivenPoint {
    driver: [GateDriverPoint; 4],
    channel_slope: [f64; 4],
}
#[derive(Clone, Copy)]
struct NodePoint {
    legs: [HalfBridgeNodePoint; 2],
    gate: Option<DrivenPoint>,
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
    fn supplied_power(self) -> Dd {
        self.gate
            .map(|g| g.driver.into_iter().fold(Dd::ZERO, |s, p| s + p.supplied))
            .unwrap_or(Dd::ZERO)
    }
    fn flow(self, out: &mut [f64]) -> Option<f64> {
        if out.len() != if self.gate.is_some() { 9 } else { 5 } {
            return None;
        }
        out.fill(0.0);
        for k in 0..2 {
            out[1 + k] = number(self.legs[k].high_current).ok()?;
            out[3 + k] = number(self.legs[k].low_current - self.legs[k].high_current).ok()?;
        }
        let mut heat = self.legs[0].loss + self.legs[1].loss;
        if let Some(g) = self.gate {
            for k in 0..4 {
                out[5 + k] = -g.driver[k].current_a;
                heat = heat + g.driver[k].loss;
            }
        }
        number(heat).ok()
    }
    fn tangent(self, de: &[f64], out: &mut [f64]) -> bool {
        let n = if self.gate.is_some() { 9 } else { 5 };
        if de.len() != n || out.len() != n {
            return false;
        }
        out.fill(0.0);
        for k in 0..2 {
            let mut high = Dd::from_f64(self.legs[k].d_high_current_d_drain_source_voltage_s)
                * (Dd::from_f64(de[1 + k]) - Dd::from_f64(de[3 + k]));
            let mut low = Dd::from_f64(self.legs[k].d_low_current_d_drain_source_voltage_s)
                * Dd::from_f64(de[3 + k]);
            if let Some(g) = self.gate {
                high = high + Dd::from_f64(g.channel_slope[2 * k]) * Dd::from_f64(de[5 + 2 * k]);
                low = low + Dd::from_f64(g.channel_slope[2 * k + 1]) * Dd::from_f64(de[6 + 2 * k]);
            }
            out[1 + k] = high.to_f64();
            out[3 + k] = (low - high).to_f64();
        }
        if let Some(g) = self.gate {
            for k in 0..4 {
                out[5 + k] =
                    (Dd::from_f64(g.driver[k].conductance_s) * Dd::from_f64(de[5 + k])).to_f64();
            }
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
            gate: None,
            legs: [
                self.legs[0].evaluate_at_node(e[1], e[3])?,
                self.legs[1].evaluate_at_node(e[2], e[4])?,
            ],
        })
    }
}

struct ConverterAccumulator {
    heat: [Dd; 13],
    work: [Dd; 2],
    absolute_work: Dd,
    parts: [Dd; 15],
    absolute_parts: [Dd; 15],
    delta: Dd,
    means: [Dd; 5],
    last: Option<([f64; 5], [DeviceConductionPoint; 4])>,
    gate_work: [Dd; 4],
    gate_voltage: [Dd; 4],
    gate_current: [Dd; 4],
    last_gate: Option<([f64; 4], [f64; 4], [f64; 4])>,
}
impl ConverterAccumulator {
    fn new() -> Self {
        Self {
            heat: [Dd::ZERO; 13],
            work: [Dd::ZERO; 2],
            absolute_work: Dd::ZERO,
            parts: [Dd::ZERO; 15],
            absolute_parts: [Dd::ZERO; 15],
            delta: Dd::ZERO,
            means: [Dd::ZERO; 5],
            last: None,
            gate_work: [Dd::ZERO; 4],
            gate_voltage: [Dd::ZERO; 4],
            gate_current: [Dd::ZERO; 4],
            last_gate: None,
        }
    }
    fn push(&mut self, a: &CapacitiveConverterAccount, dt: f64) {
        for (s, v) in self.heat.iter_mut().zip(a.all_heat_joules()) {
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
        if let Some(g) = &a.gate {
            for k in 0..4 {
                self.gate_work[k] = self.gate_work[k] + Dd::from_f64(g.source_work_j[k]);
                self.absolute_work = self.absolute_work + Dd::from_f64(g.source_work_j[k].abs());
                self.gate_voltage[k] =
                    self.gate_voltage[k] + Dd::from_f64(dt) * Dd::from_f64(g.mean_voltage_v[k]);
                self.gate_current[k] =
                    self.gate_current[k] + Dd::from_f64(dt) * Dd::from_f64(g.mean_current_a[k]);
            }
            self.last_gate = Some((g.mean_voltage_v, g.mean_current_a, g.last_driver_current_a));
        }
        self.last = Some((means, a.last_devices));
    }
    fn finish<const N: usize>(
        self,
        core: &ConverterCore<N>,
        initial: &[f64; N],
        next: &[f64; N],
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
            gate: self
                .last_gate
                .map(|(voltage, current, last)| GatePowerAccount {
                    gate_source_energy_change_j: core::array::from_fn(|k| s[7 + k]),
                    gate_drain_energy_change_j: core::array::from_fn(|k| s[11 + k]),
                    source_work_j: self.gate_work.map(Dd::to_f64),
                    resistor_heat_j: [heat[9], heat[10], heat[11], heat[12]],
                    mean_voltage_v: if plan.integration_steps == 1 {
                        voltage
                    } else {
                        self.gate_voltage
                            .map(|v| (v / Dd::from_f64(plan.period_s)).to_f64())
                    },
                    mean_current_a: if plan.integration_steps == 1 {
                        current
                    } else {
                        self.gate_current
                            .map(|v| (v / Dd::from_f64(plan.period_s)).to_f64())
                    },
                    last_driver_current_a: last,
                }),
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
        for k in 0..15 {
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
