//! Static device commutation with explicitly authored leakage and body diodes.
//!
//! Four independently commanded, finite-resistance MOS channels and four
//! unilateral finite-slope body diodes close both half-bridge midpoint KCLs.
//! A strictly positive OFF conductance is physical caller data, never an
//! invented numerical regularizer. It gives a unique algebraic node even at
//! zero current or during both-off dead time. Both-on commands describe actual
//! finite-resistance shoot-through and must remain inside authored device limits.
//!
//! The base lane retains only the original inductor flux and two bus charges.
//! [`capacitance`] adds explicit device charge and stored switch-node voltages.
//! [`gate`] supplies original gate laws and finite driver branches;
//! [`capacitance::gate`] retains reciprocal gate/Miller storage and transitions.
//! Reverse recovery and avalanche are absent. Refine steps to resolve diode
//! branches and the actual RC dynamics.
//! No default device data, efficiency factor or unobserved switching loss is used.

use super::power_load::{
    PowerLoadControls, PowerLoadPoint, RegulatedPowerLoad, VoltageSources, load_point,
};
use super::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters, invalid, poll, positive};
use crate::{NonlinearPortPower, PhsError, PreparedStepError, PreparedStepRecord, StepWorkspace};
use fs_math::dd::Dd;

/// Original finite battery and motor shafts with these same device laws.
pub mod battery_motor;

/// Finite output-capacitor charge and energy through physical commutation.
pub mod capacitance;

/// Source-owned gate-voltage channel laws and finite-impedance driver branches.
pub mod gate;

mod cycle;
pub use cycle::PreparedCommutatedSwitchCycle;

/// Centered complementary PWM with explicit non-overlap and physical gate timing.
pub mod pwm;

/// Authored static forward law I=max(0,(V_source-drain−Vf)/Rd).
/// Reverse current is carried only by the separate OFF/ON channel resistance;
/// no avalanche, recovery charge or exponential leakage is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BodyDiode {
    pub forward_voltage_v: f64,
    pub differential_resistance_ohm: f64,
}

/// Actual terminal/conduction support, enforced before and after every substep
/// and at its discrete effort. These are source limits, not current clipping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeviceOperatingBounds {
    /// Positive drain-to-source voltage, i.e. reverse bias of the body diode.
    pub maximum_reverse_voltage_v: f64,
    pub maximum_channel_current_a: f64,
    pub maximum_diode_current_a: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwitchDevice {
    pub on_resistance_ohm: f64,
    /// Finite, explicitly authored OFF leakage resistance. There is no default.
    pub off_resistance_ohm: f64,
    pub body_diode: BodyDiode,
    pub bounds: DeviceOperatingBounds,
}
impl SwitchDevice {
    pub fn validate(&self) -> Result<(), PhsError> {
        let d = self.body_diode;
        if !d.forward_voltage_v.is_finite()
            || d.forward_voltage_v < 0.0
            || [
                self.on_resistance_ohm,
                self.off_resistance_ohm,
                d.differential_resistance_ohm,
            ]
            .iter()
            .any(|r| !positive(*r) || !positive(1.0 / r))
            || [
                self.bounds.maximum_reverse_voltage_v,
                self.bounds.maximum_channel_current_a,
                self.bounds.maximum_diode_current_a,
            ]
            .iter()
            .any(|x| !positive(*x))
            || !((1.0 / self.on_resistance_ohm) + (1.0 / d.differential_resistance_ohm)).is_finite()
            || !((1.0 / self.off_resistance_ohm) + (1.0 / d.differential_resistance_ohm))
                .is_finite()
        {
            return Err(invalid(
                "commutation requires authored finite positive resistances/bounds and nonnegative diode drop",
            ));
        }
        Ok(())
    }
}

/// Device order throughout is input high, input low, output high, output low.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CommutationParameters {
    pub inductance_h: f64,
    pub bus_capacitance_f: [f64; 2],
    pub winding_resistance_ohm: f64,
    pub devices: [SwitchDevice; 4],
}
impl CommutationParameters {
    fn storage_parameters(self) -> BuckBoostParameters {
        BuckBoostParameters {
            inductance_h: self.inductance_h,
            bus_capacitance_f: self.bus_capacitance_f,
            winding_resistance_ohm: self.winding_resistance_ohm,
            input_on_resistance_ohm: [0.0; 2],
            output_on_resistance_ohm: [0.0; 2],
        }
    }
    fn winding_model(self) -> Result<AveragedBuckBoost, PhsError> {
        AveragedBuckBoost::new(
            self.storage_parameters(),
            BridgeModulation {
                input_high_fraction: 0.0,
                output_high_fraction: 0.0,
            },
        )
    }
    pub fn validate(&self) -> Result<(), PhsError> {
        for device in self.devices {
            device.validate()?;
        }
        self.winding_model()?;
        Ok(())
    }
    /// Construct the SAME original [inductor flux, input charge, output charge].
    /// Device operating support additionally depends on the actual gate command.
    pub fn state(&self, current_a: f64, bus_voltage_v: [f64; 2]) -> Result<[f64; 3], PhsError> {
        self.validate()?;
        self.winding_model()?.state(current_a, bus_voltage_v)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DeviceGateState {
    /// Independent input-high, input-low, output-high, output-low commands.
    pub on: [bool; 4],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CommutationInterval {
    pub gates: DeviceGateState,
    pub duration_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeviceConductionPoint {
    pub drain_source_voltage_v: f64,
    /// Signed channel current from drain to source, including declared OFF leak.
    pub channel_current_a: f64,
    /// Nonnegative body-diode current from source to drain.
    pub diode_current_a: f64,
    pub terminal_current_a: f64,
    pub channel_power_w: f64,
    pub diode_power_w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfBridgePoint {
    pub switch_node_voltage_v: f64,
    pub high: DeviceConductionPoint,
    pub low: DeviceConductionPoint,
    pub d_node_voltage_d_bus_voltage: f64,
    /// Signed derivative with respect to current DRAWN FROM the switch node.
    pub d_node_voltage_d_draw_current_ohm: f64,
    pub d_high_current_d_bus_voltage_s: f64,
    pub d_high_current_d_draw_current: f64,
    node: Dd,
    high_current: Dd,
    loss: Dd,
}

/// Conduction at a physically stored switch-node voltage. The difference
/// `high.terminal_current_a - low.terminal_current_a` is the current delivered
/// to the node; its capacitor and connected branches own the node KCL. This
/// evaluation never replaces that voltage with an algebraic equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfBridgeNodePoint {
    pub switch_node_voltage_v: f64,
    pub high: DeviceConductionPoint,
    pub low: DeviceConductionPoint,
    /// Incremental terminal conductance of the high device at its actual Vds.
    pub d_high_current_d_drain_source_voltage_s: f64,
    /// Incremental terminal conductance of the low device at its actual Vds.
    pub d_low_current_d_drain_source_voltage_s: f64,
    high_current: Dd,
    low_current: Dd,
    loss: Dd,
}

/// Reusable bounded static KCL owner. High current enters from the bus; low
/// current leaves to the return; I_high−I_low=current drawn from the node.
#[derive(Clone, Copy, Debug)]
pub struct HalfBridgeConduction {
    devices: [SwitchDevice; 2],
    conductance_s: [f64; 2],
}
fn number(value: Dd) -> Result<f64, PhsError> {
    if !value.hi.is_finite() || !value.lo.is_finite() || !value.to_f64().is_finite() {
        Err(invalid("commutation arithmetic is unrepresentable"))
    } else {
        Ok(value.to_f64())
    }
}
fn total(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
        .to_f64()
}

impl HalfBridgeConduction {
    pub fn new(high: SwitchDevice, low: SwitchDevice, on: [bool; 2]) -> Result<Self, PhsError> {
        let devices = [high, low];
        let conductance_s = core::array::from_fn(|i| {
            1.0 / if on[i] {
                devices[i].on_resistance_ohm
            } else {
                devices[i].off_resistance_ohm
            }
        });
        Self::with_conductance(high, low, conductance_s)
    }

    /// Reuse the original physical diode law and operating support with an
    /// independently admitted actual channel conductance. Driven channels use
    /// their complete G(Vgs) law; the command bit does not select Ron/Roff.
    pub(super) fn with_conductance(
        high: SwitchDevice,
        low: SwitchDevice,
        conductance_s: [f64; 2],
    ) -> Result<Self, PhsError> {
        high.validate()?;
        low.validate()?;
        let devices = [high, low];
        if conductance_s.iter().enumerate().any(|(i, &g)| {
            !positive(g)
                || !(g + 1.0 / devices[i].body_diode.differential_resistance_ohm).is_finite()
        }) {
            return Err(invalid(
                "positive representable original channel conductances required",
            ));
        }
        Ok(Self {
            devices,
            conductance_s,
        })
    }

    fn device_point(
        &self,
        index: usize,
        voltage: Dd,
    ) -> Result<(DeviceConductionPoint, Dd, f64, Dd), PhsError> {
        let v = number(voltage)?;
        let channel = voltage * Dd::from_f64(self.conductance_s[index]);
        let channel_current_a = number(channel)?;
        if v != 0.0 && channel_current_a == 0.0 {
            return Err(invalid("nonzero channel current is unrepresentable"));
        }
        self.device_point_with_channel(index, voltage, channel, self.conductance_s[index])
    }

    /// One diode/terminal/loss owner serves both the original linear channel and
    /// an actual nonlinear channel current. The latter enters as its represented
    /// state-flow current, so its power uses exactly the same physical current.
    fn device_point_with_channel(
        &self,
        index: usize,
        voltage: Dd,
        channel: Dd,
        channel_derivative_s: f64,
    ) -> Result<(DeviceConductionPoint, Dd, f64, Dd), PhsError> {
        let device = self.devices[index];
        let v = number(voltage)?;
        let channel_current_a = number(channel)?;
        if !channel_derivative_s.is_finite() {
            return Err(invalid("finite original channel drain derivative required"));
        }
        let overdrive = -voltage - Dd::from_f64(device.body_diode.forward_voltage_v);
        let forward = number(overdrive)? > 0.0;
        let diode = if forward {
            overdrive / Dd::from_f64(device.body_diode.differential_resistance_ohm)
        } else {
            Dd::ZERO
        };
        let diode_current_a = number(diode)?;
        if forward && diode_current_a <= 0.0 {
            return Err(invalid("nonzero body-diode current is unrepresentable"));
        }
        let channel_power = voltage * channel;
        let diode_power = -voltage * diode;
        let channel_power_w = number(channel_power)?;
        let diode_power_w = number(diode_power)?;
        if channel_power_w < 0.0
            || diode_power_w < 0.0
            || (channel_current_a != 0.0 && channel_power_w == 0.0)
            || (diode_current_a != 0.0 && diode_power_w == 0.0)
        {
            return Err(invalid("positive physical device power is unrepresentable"));
        }
        let current = channel - diode;
        let slope = number(
            Dd::from_f64(channel_derivative_s)
                + Dd::from_f64(if forward {
                    1.0 / device.body_diode.differential_resistance_ohm
                } else {
                    0.0
                }),
        )?;
        Ok((
            DeviceConductionPoint {
                drain_source_voltage_v: v,
                channel_current_a,
                diode_current_a,
                terminal_current_a: number(current)?,
                channel_power_w,
                diode_power_w,
            },
            current,
            slope,
            channel_power + diode_power,
        ))
    }

    fn support(&self, index: usize, point: DeviceConductionPoint) -> Result<(), PhsError> {
        let b = self.devices[index].bounds;
        if point.drain_source_voltage_v > b.maximum_reverse_voltage_v
            || point.channel_current_a.abs() > b.maximum_channel_current_a
            || point.diode_current_a > b.maximum_diode_current_a
        {
            return Err(invalid(
                "actual commutated device voltage/current leaves authored support",
            ));
        }
        Ok(())
    }

    /// Evaluate the original channel and diode laws at an independent stored
    /// node voltage. Actual voltage/current support is enforced for both
    /// devices, including both-off and both-on commands. No current clipping,
    /// node reset, capacitive loss or energy debit is introduced here.
    pub fn evaluate_at_node(
        &self,
        bus_voltage_v: f64,
        switch_node_voltage_v: f64,
    ) -> Result<HalfBridgeNodePoint, PhsError> {
        if !bus_voltage_v.is_finite() || !switch_node_voltage_v.is_finite() {
            return Err(invalid(
                "finite half-bridge bus and stored node voltages required",
            ));
        }
        self.evaluate_at_node_dd(
            Dd::from_f64(bus_voltage_v),
            Dd::from_f64(switch_node_voltage_v),
        )
    }

    /// Evaluate two original channel-only current responses at a stored node.
    /// Diode response, terminal support and physical heat are the same owners
    /// used by prescribed-conductance commutation. A saturation plateau may have
    /// zero drain derivative; no channel current is divided by drain voltage.
    pub(super) fn evaluate_at_node_with_channels(
        &self,
        bus_voltage_v: f64,
        switch_node_voltage_v: f64,
        channels: [gate::GateChannelPoint; 2],
    ) -> Result<HalfBridgeNodePoint, PhsError> {
        if !bus_voltage_v.is_finite() || !switch_node_voltage_v.is_finite() {
            return Err(invalid(
                "finite half-bridge bus and stored node voltages required",
            ));
        }
        let node = Dd::from_f64(switch_node_voltage_v);
        let high_voltage = Dd::from_f64(bus_voltage_v) - node;
        let high = channels[0].admit(number(high_voltage)?)?;
        let low = channels[1].admit(switch_node_voltage_v)?;
        self.node_point_from_devices(
            switch_node_voltage_v,
            self.device_point_with_channel(
                0,
                high_voltage,
                Dd::from_f64(high.current_a),
                high.drain_derivative_s,
            )?,
            self.device_point_with_channel(
                1,
                node,
                Dd::from_f64(low.current_a),
                low.drain_derivative_s,
            )?,
        )
    }

    fn evaluate_at_node_dd(
        &self,
        bus_voltage_v: Dd,
        switch_node_voltage_v: Dd,
    ) -> Result<HalfBridgeNodePoint, PhsError> {
        number(bus_voltage_v)?;
        let node = number(switch_node_voltage_v)?;
        self.node_point_from_devices(
            node,
            self.device_point(0, bus_voltage_v - switch_node_voltage_v)?,
            self.device_point(1, switch_node_voltage_v)?,
        )
    }

    fn node_point_from_devices(
        &self,
        node: f64,
        high: (DeviceConductionPoint, Dd, f64, Dd),
        low: (DeviceConductionPoint, Dd, f64, Dd),
    ) -> Result<HalfBridgeNodePoint, PhsError> {
        let (high, high_current, high_slope, high_loss) = high;
        let (low, low_current, low_slope, low_loss) = low;
        self.support(0, high)?;
        self.support(1, low)?;
        Ok(HalfBridgeNodePoint {
            switch_node_voltage_v: node,
            high,
            low,
            d_high_current_d_drain_source_voltage_s: high_slope,
            d_low_current_d_drain_source_voltage_s: low_slope,
            high_current,
            low_current,
            loss: high_loss + low_loss,
        })
    }

    /// Enumerate four affine diode active sets and admit only a solution of the
    /// ORIGINAL unilateral laws and KCL within arithmetic roundoff. This is a
    /// fixed amount of work, not a second iterative dynamics solver. Negative
    /// bus voltage can forward-bias both diodes and is retained when supported.
    pub fn evaluate(
        &self,
        bus_voltage_v: f64,
        draw_current_a: f64,
    ) -> Result<HalfBridgePoint, PhsError> {
        if !bus_voltage_v.is_finite() || !draw_current_a.is_finite() {
            return Err(invalid(
                "finite half-bridge bus voltage and drawn current required",
            ));
        }
        let g = self.conductance_s;
        let gd = self
            .devices
            .map(|d| 1.0 / d.body_diode.differential_resistance_ohm);
        let vf = self.devices.map(|d| d.body_diode.forward_voltage_v);
        for active in 0..4 {
            let high_on = active & 1 != 0;
            let low_on = active & 2 != 0;
            let a = g[0] + if high_on { gd[0] } else { 0.0 };
            let b = g[1] + if low_on { gd[1] } else { 0.0 };
            let scale = a.max(b);
            let an = a / scale;
            let bn = b / scale;
            let numerator = Dd::from_f64(an) * Dd::from_f64(bus_voltage_v)
                + Dd::from_f64(if high_on { gd[0] / scale } else { 0.0 }) * Dd::from_f64(vf[0])
                - Dd::from_f64(if low_on { gd[1] / scale } else { 0.0 }) * Dd::from_f64(vf[1])
                - Dd::from_f64(draw_current_a) / Dd::from_f64(scale);
            let node = numerator / (Dd::from_f64(an) + Dd::from_f64(bn));
            let Ok((high, high_current, ah, high_loss)) =
                self.device_point(0, Dd::from_f64(bus_voltage_v) - node)
            else {
                continue;
            };
            let Ok((low, low_current, bl, low_loss)) = self.device_point(1, node) else {
                continue;
            };
            let defect = number(high_current - low_current - Dd::from_f64(draw_current_a))?;
            let magnitude = total([
                high.terminal_current_a.abs(),
                low.terminal_current_a.abs(),
                draw_current_a.abs(),
            ]);
            let tolerance = 64.0 * f64::EPSILON * magnitude;
            if !magnitude.is_finite() || !tolerance.is_finite() || defect.abs() > tolerance {
                continue;
            }
            self.support(0, high)?;
            self.support(1, low)?;
            let largest = ah.max(bl);
            let denominator = ah / largest + bl / largest;
            let fraction = (ah / largest) / denominator;
            let resistance = (1.0 / largest) / denominator;
            let smallest = ah.min(bl);
            let parallel = smallest / (1.0 + smallest / largest);
            if !positive(fraction) || !positive(resistance) || !positive(parallel) {
                return Err(invalid(
                    "half-bridge exact differential response is unrepresentable",
                ));
            }
            return Ok(HalfBridgePoint {
                switch_node_voltage_v: number(node)?,
                high,
                low,
                d_node_voltage_d_bus_voltage: fraction,
                d_node_voltage_d_draw_current_ohm: -resistance,
                d_high_current_d_bus_voltage_s: parallel,
                d_high_current_d_draw_current: fraction,
                node,
                high_current,
                loss: high_loss + low_loss,
            });
        }
        Err(invalid(
            "no representable half-bridge solution satisfies original diode laws and KCL",
        ))
    }
}

/// Physical native receipt, distinct from the ideal converter's five-only heat
/// record. Native dissipation contains winding, four channels and four diodes;
/// source-resistor heat and regulated load work retain separate ownership.
#[derive(Clone, Debug)]
pub struct CommutatedStep {
    pub native: PreparedStepRecord,
    pub terminal_work_j: [f64; 2],
    pub winding_heat_j: f64,
    pub channel_heat_j: [f64; 4],
    pub diode_heat_j: [f64; 4],
    pub source_work_j: [f64; 2],
    pub terminal_resistor_heat_j: [f64; 2],
    pub load_work_j: [f64; 2],
    pub curtailed_load_energy_j: [f64; 2],
    pub mean_source_current_a: [f64; 2],
    /// Actual final discrete effort only; endpoint/substep limits are checked
    /// internally throughout the cycle, never on this last sample alone.
    pub last_devices: [DeviceConductionPoint; 4],
    pub last_loads: [PowerLoadPoint; 2],
    pub integration_steps: usize,
}

struct OperatingPoint {
    legs: [HalfBridgePoint; 2],
    loads: [PowerLoadPoint; 2],
    source_current_a: [f64; 2],
}
struct PreparedInterval {
    model: AveragedBuckBoost,
    legs: [HalfBridgeConduction; 2],
    workspace: StepWorkspace,
    inverse_storage: [f64; 3],
    duration_s: f64,
    controls: PowerLoadControls,
}
impl PreparedInterval {
    fn new(
        parameters: CommutationParameters,
        gates: DeviceGateState,
        duration_s: f64,
        controls: PowerLoadControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(cancelled)?;
        parameters.validate()?;
        if !positive(duration_s)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid(
                "commutation requires positive duration and finite native energy controls",
            )
            .into());
        }
        controls.solve.validate()?;
        let legs = [
            HalfBridgeConduction::new(
                parameters.devices[0],
                parameters.devices[1],
                [gates.on[0], gates.on[1]],
            )?,
            HalfBridgeConduction::new(
                parameters.devices[2],
                parameters.devices[3],
                [gates.on[2], gates.on[3]],
            )?,
        ];
        let model = parameters.winding_model()?;
        let mut workspace = StepWorkspace::new(model.system())?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        let inverse_storage = [
            1.0 / parameters.inductance_h,
            1.0 / parameters.bus_capacitance_f[0],
            1.0 / parameters.bus_capacitance_f[1],
        ];
        poll(cancelled)?;
        Ok(Self {
            model,
            legs,
            workspace,
            inverse_storage,
            duration_s,
            controls,
        })
    }

    fn point(
        legs: &[HalfBridgeConduction; 2],
        effort: &[f64],
        sources: VoltageSources,
        loads: &[Option<RegulatedPowerLoad>; 2],
    ) -> Result<OperatingPoint, PhsError> {
        if effort.len() != 3 || effort.iter().any(|e| !e.is_finite()) {
            return Err(invalid("finite original commutation efforts required"));
        }
        let legs = [
            legs[0].evaluate(effort[1], effort[0])?,
            legs[1].evaluate(effort[2], -effort[0])?,
        ];
        let loads = [
            load_point(loads[0], effort[1])?,
            load_point(loads[1], effort[2])?,
        ];
        let source_current_a = core::array::from_fn(|p| sources.current_a(p, effort[p + 1]));
        if source_current_a.iter().any(|i| !i.is_finite()) {
            return Err(invalid("commutation source current is unrepresentable"));
        }
        Ok(OperatingPoint {
            legs,
            loads,
            source_current_a,
        })
    }

    fn step(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[f64; 2],
        terminal_resistance_ohm: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; 3],
        voltage: &mut [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CommutatedStep, PreparedStepError> {
        poll(cancelled)?;
        let sources = VoltageSources::new(*source_voltage_v, *terminal_resistance_ohm)?;
        let initial_effort: [f64; 3] =
            core::array::from_fn(|i| initial[i] * self.inverse_storage[i]);
        Self::point(&self.legs, &initial_effort, sources, loads)?;
        let legs = self.legs;
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if out.len() != 3 {
                return None;
            }
            let point = Self::point(&legs, effort, sources, loads).ok()?;
            out[0] = number(point.legs[1].node - point.legs[0].node).ok()?;
            let mut supplied = Dd::ZERO;
            for p in 0..2 {
                let external = Dd::from_f64(point.source_current_a[p])
                    - Dd::from_f64(point.loads[p].current_a);
                out[p + 1] = number(point.legs[p].high_current - external).ok()?;
                supplied = supplied + Dd::from_f64(effort[p + 1]) * external;
            }
            Some(NonlinearPortPower {
                supplied_power_w: number(supplied).ok()?,
                dissipated_power_w: number(point.legs[0].loss + point.legs[1].loss).ok()?,
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
            if de.len() != 3 || out.len() != 3 {
                return false;
            }
            let Ok(point) = Self::point(&legs, effort, sources, loads) else {
                return false;
            };
            let mut node = [0.0; 2];
            for p in 0..2 {
                let ds = if p == 0 { de[0] } else { -de[0] };
                let leg = point.legs[p];
                node[p] = total([
                    leg.d_node_voltage_d_bus_voltage * de[p + 1],
                    leg.d_node_voltage_d_draw_current_ohm * ds,
                ]);
                out[p + 1] = total([
                    leg.d_high_current_d_bus_voltage_s * de[p + 1],
                    leg.d_high_current_d_draw_current * ds,
                    (point.loads[p].incremental_conductance_s + sources.conductance_s[p])
                        * de[p + 1],
                ]);
            }
            out[0] = total([node[1], -node[0]]);
            out.iter().all(|v| v.is_finite())
        };
        let (mut candidate, mut effort_v) = ([0.0; 3], [0.0; 2]);
        let native = self.workspace.step_into_supplied_controlled(
            self.model.system(),
            initial,
            &[0.0; 2],
            self.duration_s,
            &mut candidate,
            &mut effort_v,
            &port,
            Some(&hessian),
            Some(&tangent),
            &mut *cancelled,
        )?;
        let effort = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("commutation requires the accepted original discrete effort"))?;
        let point = Self::point(&self.legs, effort, sources, loads)?;
        let final_effort: [f64; 3] =
            core::array::from_fn(|i| candidate[i] * self.inverse_storage[i]);
        Self::point(&self.legs, &final_effort, sources, loads)?;
        let devices = [
            point.legs[0].high,
            point.legs[0].low,
            point.legs[1].high,
            point.legs[1].low,
        ];
        let h = Dd::from_f64(self.duration_s);
        let winding_heat_j = number(
            h * Dd::from_f64(self.model.parameters.winding_resistance_ohm)
                * Dd::from_f64(effort[0])
                * Dd::from_f64(effort[0]),
        )?;
        let channel_heat_j = devices.map(|d| (h * Dd::from_f64(d.channel_power_w)).to_f64());
        let diode_heat_j = devices.map(|d| (h * Dd::from_f64(d.diode_power_w)).to_f64());
        let load_work_j = point
            .loads
            .map(|load| (h * Dd::from_f64(load.delivered_power_w)).to_f64());
        let curtailed_load_energy_j = point
            .loads
            .map(|load| (h * Dd::from_f64(load.curtailed_power_w)).to_f64());
        let source_work_j = core::array::from_fn(|p| {
            (h * Dd::from_f64(source_voltage_v[p]) * Dd::from_f64(point.source_current_a[p]))
                .to_f64()
        });
        let terminal_resistor_heat_j = core::array::from_fn(|p| {
            (h * Dd::from_f64(terminal_resistance_ohm[p])
                * Dd::from_f64(point.source_current_a[p])
                * Dd::from_f64(point.source_current_a[p]))
            .to_f64()
        });
        let terminal_work_j = core::array::from_fn(|p| {
            (h * Dd::from_f64(effort_v[p]) * Dd::from_f64(point.source_current_a[p])
                - Dd::from_f64(load_work_j[p]))
            .to_f64()
        });
        let heat = total(
            [winding_heat_j]
                .into_iter()
                .chain(channel_heat_j)
                .chain(diode_heat_j),
        );
        let scale = total(source_work_j.map(f64::abs))
            .max(total(load_work_j))
            .max(total(terminal_resistor_heat_j))
            .max(native.delta_h.abs())
            .max(heat);
        let allowance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            total([native.delta_h, heat, -native.supplied]),
            total([native.dissipated, -heat]),
            total(terminal_work_j.into_iter().chain([-native.supplied])),
            total(
                [native.delta_h, heat]
                    .into_iter()
                    .chain(terminal_resistor_heat_j)
                    .chain(load_work_j)
                    .chain(source_work_j.map(|x| -x)),
            ),
        ];
        if !scale.is_finite()
            || !allowance.is_finite()
            || !winding_heat_j.is_finite()
            || winding_heat_j < 0.0
            || channel_heat_j
                .iter()
                .chain(&diode_heat_j)
                .chain(&terminal_resistor_heat_j)
                .chain(&load_work_j)
                .chain(&curtailed_load_energy_j)
                .any(|x| !x.is_finite() || *x < 0.0)
            || source_work_j
                .iter()
                .chain(&terminal_work_j)
                .any(|x| !x.is_finite())
            || residuals
                .iter()
                .any(|r| !r.is_finite() || r.abs() > allowance)
        {
            return Err(invalid(
                "original commutation storage, device heat and source/load work do not close",
            )
            .into());
        }
        poll(cancelled)?;
        next.copy_from_slice(&candidate);
        voltage.copy_from_slice(&effort_v);
        Ok(CommutatedStep {
            native,
            terminal_work_j,
            winding_heat_j,
            channel_heat_j,
            diode_heat_j,
            source_work_j,
            terminal_resistor_heat_j,
            load_work_j,
            curtailed_load_energy_j,
            mean_source_current_a: point.source_current_a,
            last_devices: devices,
            last_loads: point.loads,
            integration_steps: 1,
        })
    }
}
