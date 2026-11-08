//! Finite battery, physical leads and a changing-duty synchronous motor drive.
//!
//! All nine original charges, winding fluxes, angle and shaft momenta evolve in
//! one native solve. The two DC capacitors retain their separate voltages and
//! inrush energy. Their current-input/voltage-output ports are joined by the
//! physical lead law `(V_battery - V_dc)/R_lead`, not a skew connection between
//! two voltage outputs. Regeneration follows the same signed current path.
//!
//! Original SOC/OCV and saturated magnetic storage remain the constitutive
//! owners. Changing duty changes only reciprocal coupling and channel loss.
//! This inherits the averaged inverter and reversible equivalent-circuit battery
//! bounds: no PWM ripple, switching heat, battery protection or thermal drift.
//! SOC exhaustion, overcharge and negative inverter bus voltage refuse without
//! clipping; endpoint admission does not certify between-endpoint bounds.

use super::{AveragedSynchronousInverter, BridgeStamp, InverterInterval, InverterParameters};
use crate::circuit::battery::{BatteryInterval, BatteryOcvCurve, BatteryParameters, BatterySource};
use crate::circuit::synchronous::{
    SynchronousMotorParameters, SynchronousOperatingPoint, invalid, poll, reference, storage,
};
use crate::{
    InterconnectionBudget, PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord,
    ResistanceWitness, SolveBudget, StepWorkspace, interconnect_mapped, storage_energy_change,
};
use fs_math::{dd::Dd, det};

pub mod shaft;

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}

/// Native nonlinear solve limits and independent physical work allowances.
/// Relative allowances scale with interval work/change, never background energy.
#[derive(Clone, Copy, Debug)]
pub struct BatteryInverterControls {
    pub solve: SolveBudget,
    pub backtracking_halvings: u8,
    pub absolute_energy_tolerance_j: f64,
    pub relative_energy_tolerance: f64,
}

/// Complete nine-state native solve and independently owned physical accounts.
/// `inverter.source_work_j` enters from the lead; it is INTERNAL to this drive.
/// Only the two shaft works enter the complete system from outside.
#[derive(Clone, Debug)]
pub struct BatteryInverterStep {
    pub native: PreparedStepRecord,
    pub battery: BatteryInterval,
    pub inverter: InverterInterval,
    pub inverter_storage_change_j: f64,
    /// Positive battery-to-inverter current; negative during regeneration [A].
    pub lead_current_a: f64,
    pub lead_heat_j: f64,
    /// Work-conjugate discrete battery/inverter bus voltages, not exact means.
    pub discrete_bus_voltage_v: [f64; 2],
}

impl BatteryInverterStep {
    /// Battery series/polarization, lead, then the inverter's eight heat owners.
    pub fn heat_joules(&self) -> [f64; 11] {
        let mut heat = [0.0; 11];
        heat[..3].copy_from_slice(&[
            self.battery.series_heat_j,
            self.battery.polarization_heat_j,
            self.lead_heat_j,
        ]);
        heat[3..].copy_from_slice(&self.inverter.heat_joules());
        heat
    }

    pub fn balance_residual_j(&self) -> f64 {
        sum(self.heat_joules().into_iter().chain([
            self.battery.storage_change_j(),
            self.inverter_storage_change_j,
            -self.inverter.shaft_work_j[0],
            -self.inverter.shaft_work_j[1],
        ]))
    }
}

#[derive(Clone, Copy)]
struct DriveStamp {
    bridge: BridgeStamp,
    resistance: [f64; 81],
    root: [f64; 81],
}

impl DriveStamp {
    fn new(battery: &BatterySource, bridge: BridgeStamp, lead_r: f64) -> Result<Self, PhsError> {
        let mut resistance = [0.0; 81];
        for i in 0..3 {
            for j in 0..3 {
                resistance[i * 9 + j] = battery.system().r[i * 3 + j];
            }
        }
        for i in 0..6 {
            for j in 0..6 {
                resistance[(i + 3) * 9 + j + 3] = bridge.resistance[i * 6 + j];
            }
        }
        let g = lead_r.recip();
        resistance[2 * 9 + 2] += g;
        resistance[3 * 9 + 3] += g;
        resistance[2 * 9 + 3] = -g;
        resistance[3 * 9 + 2] = -g;
        // Retain every nonzero original source row. The battery has two rows,
        // the inverter at most three and the lead one; no QR or new factorization
        // is needed when modulation changes. Zero rows carry no physical loss.
        let mut root = [0.0; 81];
        let mut row = 0;
        let battery_rows = battery
            .system()
            .resistance_witness()
            .ok_or_else(|| invalid("finite battery must retain its source loss witness"))?
            .rows();
        for (rows, width, offset) in [(battery_rows, 3, 0), (bridge.root.as_slice(), 6, 3)] {
            for source in rows.chunks_exact(width) {
                if source.iter().any(|v| *v != 0.0) {
                    if row >= 8 {
                        return Err(invalid("battery drive source rows exceed native dimension"));
                    }
                    root[row * 9 + offset..row * 9 + offset + width].copy_from_slice(source);
                    row += 1;
                }
            }
        }
        root[row * 9 + 2] = det::sqrt(g);
        root[row * 9 + 3] = -det::sqrt(g);
        if resistance.iter().chain(&root).any(|v| !v.is_finite()) {
            return Err(invalid(
                "battery inverter resistance or source rows are unrepresentable",
            ));
        }
        crate::linear_exact::positive_heat::check_resistance_root(&resistance, &root, 9)?;
        Ok(Self {
            bridge,
            resistance,
            root,
        })
    }
}

/// State order: battery chemical/polarization/bus charges, inverter bus charge,
/// stator alpha/beta fluxes, relative shaft angle, rotor and stator momenta.
pub struct BatterySynchronousInverter {
    battery: BatterySource,
    inverter: AveragedSynchronousInverter,
    lead_resistance_ohm: f64,
    system: PortHamiltonian,
    held: DriveStamp,
}

impl AveragedSynchronousInverter {
    /// Join both original finite buses through a positive physical lead. Zero
    /// resistance requires a common-voltage constraint and is not approximated.
    pub fn with_battery_source(
        mut self,
        mut battery: BatterySource,
        lead_resistance_ohm: f64,
    ) -> Result<BatterySynchronousInverter, PhsError> {
        if !lead_resistance_ohm.is_finite()
            || lead_resistance_ohm <= 0.0
            || !lead_resistance_ohm.recip().is_finite()
            || lead_resistance_ohm.recip() <= 0.0
        {
            return Err(invalid(
                "battery-inverter lead resistance must be positive and representable",
            ));
        }
        let held = DriveStamp::new(&battery, self.held, lead_resistance_ohm)?;
        let mut system = interconnect_mapped(
            battery.share_system(),
            self.system.share_storage(),
            &[],
            &[],
            &[],
            InterconnectionBudget {
                max_states: 9,
                max_total_ports: 4,
                max_dense_bytes: 1 << 18,
                max_coupling_multiply_terms: 256,
            },
        )?;
        system.r.copy_from_slice(&held.resistance);
        system = system.with_resistance_witness(ResistanceWitness::try_from_source_rows(
            9,
            held.resistance.to_vec(),
            held.root.to_vec(),
        )?)?;
        Ok(BatterySynchronousInverter {
            battery,
            inverter: self,
            lead_resistance_ohm,
            system,
            held,
        })
    }
}

impl BatterySynchronousInverter {
    /// Original additive ports are battery bus current, inverter bus current,
    /// rotor torque and stator torque. The prepared owner holds both currents at
    /// zero; battery/lead exchange is internal to the original resistance law.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn parameters(&self) -> InverterParameters {
        self.inverter.parameters()
    }
    pub fn machine_parameters(&self) -> SynchronousMotorParameters {
        self.inverter.machine_parameters()
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.battery.parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.battery.curve()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }

    /// Initialize both finite bus charges independently; a mismatch drives real
    /// lead inrush. `machine_state` retains the complete original five states.
    pub fn state(
        &self,
        soc: f64,
        polarization_v: f64,
        battery_bus_v: f64,
        inverter_bus_v: f64,
        machine_state: &[f64; 5],
    ) -> Result<[f64; 9], PhsError> {
        let battery = self.battery.state(soc, polarization_v, battery_bus_v)?;
        let inverter = self.inverter.state(inverter_bus_v, machine_state)?;
        let mut state = [0.0; 9];
        state[..3].copy_from_slice(&battery);
        state[3..].copy_from_slice(&inverter);
        self.admit(&state)?;
        Ok(state)
    }

    fn admit(&self, state: &[f64; 9]) -> Result<(), PhsError> {
        self.battery.curve().voltage_v(state[0])?;
        self.inverter
            .admit_bus(&[state[3], state[4], state[5], state[6], state[7], state[8]])?;
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(state).is_finite() {
            return Err(invalid(
                "battery inverter requires complete supported finite storage",
            ));
        }
        Ok(())
    }

    fn apply_stamp(&mut self, stamp: &DriveStamp) {
        self.inverter.apply_stamp(&stamp.bridge);
        for i in 0..6 {
            for j in 0..6 {
                self.system.j[(i + 3) * 9 + j + 3] = self.inverter.system.j[i * 6 + j];
            }
        }
        self.system.r.copy_from_slice(&stamp.resistance);
        let witness = self
            .system
            .resistance_witness
            .as_mut()
            .expect("battery inverter retains its original source witness");
        witness.matrix.copy_from_slice(&stamp.resistance);
        witness.root.copy_from_slice(&stamp.root);
    }

    pub fn prepare_controlled(
        self,
        controls: BatteryInverterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedBatterySynchronousInverter, PreparedStepError> {
        poll(&mut cancelled)?;
        if !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(
                invalid("battery inverter requires finite nonnegative energy allowances").into(),
            );
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedBatterySynchronousInverter {
            model: self,
            workspace,
            controls,
        })
    }
}

/// Reusable nine-state native scratch. Warm duty changes rebuild neither the
/// interconnection nor storage and allocate nothing on successful steps.
pub struct PreparedBatterySynchronousInverter {
    model: BatterySynchronousInverter,
    workspace: StepWorkspace,
    controls: BatteryInverterControls,
}

impl PreparedBatterySynchronousInverter {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn parameters(&self) -> InverterParameters {
        self.model.parameters()
    }
    pub fn machine_parameters(&self) -> SynchronousMotorParameters {
        self.model.machine_parameters()
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.model.battery_parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.model.battery_curve()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.model.lead_resistance_ohm()
    }
    pub fn high_fraction(&self) -> [f64; 3] {
        self.model.held.bridge.high_fraction
    }

    pub fn machine_operating_point(
        &self,
        state: &[f64; 9],
    ) -> Result<SynchronousOperatingPoint, PhsError> {
        self.model.admit(state)?;
        storage::operating_point(
            self.model.inverter.machine_parameters,
            &self.model.inverter.machine_magnetic_law,
            &[state[4], state[5], state[6], state[7], state[8]],
        )
    }

    pub fn current_reference_for_torque_controlled(
        &self,
        requested_torque_nm: f64,
        d_current_a: f64,
        current_limit_a: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<reference::TorqueCurrentReference, PreparedStepError> {
        reference::allocate(
            self.model.inverter.machine_parameters,
            &self.model.inverter.machine_magnetic_law,
            requested_torque_nm,
            d_current_a,
            current_limit_a,
            cancelled,
        )
    }

    /// Hold physical duty and external shaft torques over this native interval.
    /// Signed source, lead, bridge and shaft work are independently checked.
    /// Any refusal restores the previous J/R/source witness and leaves `next`
    /// unchanged; accepted state is never projected at a duty command edge.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 9],
        high_fraction: [f64; 3],
        shaft_torque_nm: [f64; 2],
        duration_s: f64,
        next: &mut [f64; 9],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryInverterStep, PreparedStepError> {
        poll(&mut cancelled)?;
        self.model.admit(initial)?;
        let inverter = &self.model.inverter;
        let bridge = BridgeStamp::new(
            inverter.parameters,
            inverter.machine_parameters,
            &inverter.base_resistance,
            inverter.scales,
            high_fraction,
        )?;
        let stamp = DriveStamp::new(&self.model.battery, bridge, self.model.lead_resistance_ohm)?;
        poll(&mut cancelled)?;
        self.model.apply_stamp(&stamp);
        let (mut candidate, mut output) = ([0.0; 9], [0.0; 4]);
        let result: Result<BatteryInverterStep, PreparedStepError> = (|| {
            let native = self.workspace.step_into_controlled(
                &self.model.system,
                initial,
                &[0.0, 0.0, shaft_torque_nm[0], shaft_torque_nm[1]],
                duration_s,
                &mut candidate,
                &mut output,
                &mut cancelled,
            )?;
            self.model.admit(&candidate)?;
            let e = self.workspace.last_step_discrete_effort().ok_or_else(|| {
                invalid("battery inverter requires accepted constitutive efforts")
            })?;
            let lead_current_a = ((Dd::from_f64(e[2]) - Dd::from_f64(e[3]))
                / Dd::from_f64(self.model.lead_resistance_ohm))
            .to_f64();
            let battery = self.model.battery.interval(
                &[initial[0], initial[1], initial[2]],
                &[candidate[0], candidate[1], candidate[2]],
                -lead_current_a,
                duration_s,
            )?;
            let inverter = self.model.inverter.interval(
                &stamp.bridge,
                &[e[3], e[4], e[5], e[6], e[7], e[8]],
                &[lead_current_a, shaft_torque_nm[0], shaft_torque_nm[1]],
                duration_s,
            )?;
            let inverter_storage_change_j = storage_energy_change(
                self.model.inverter.system.storage(),
                &initial[3..],
                &candidate[3..],
                &mut [0.0; 6],
                &mut [0.0; 6],
            );
            let lead_heat_j =
                duration_s * lead_current_a * self.model.lead_resistance_ohm * lead_current_a;
            let step = BatteryInverterStep {
                native,
                battery,
                inverter,
                inverter_storage_change_j,
                lead_current_a,
                lead_heat_j,
                discrete_bus_voltage_v: [e[2], e[3]],
            };
            let inverter_heat = sum(inverter.heat_joules());
            let heat = sum(step.heat_joules());
            let shaft_work = sum(inverter.shaft_work_j);
            let source_scale = sum([
                battery.chemical_energy_change_j.abs(),
                battery.polarization_energy_change_j.abs(),
                battery.bus_energy_change_j.abs(),
                battery.terminal_work_j.abs(),
                battery.heat_j(),
            ]);
            let inverter_scale = sum([
                inverter_storage_change_j.abs(),
                inverter_heat,
                inverter.source_work_j.abs(),
                inverter.shaft_work_j[0].abs(),
                inverter.shaft_work_j[1].abs(),
            ]);
            let lead_scale = sum([
                battery.terminal_work_j.abs(),
                inverter.source_work_j.abs(),
                lead_heat_j,
            ]);
            let balances = [
                (battery.balance_residual_j(), source_scale),
                (
                    sum([
                        inverter_storage_change_j,
                        inverter_heat,
                        -inverter.source_work_j,
                        -shaft_work,
                    ]),
                    inverter_scale,
                ),
                (
                    sum([battery.terminal_work_j, inverter.source_work_j, lead_heat_j]),
                    lead_scale,
                ),
                (
                    sum([
                        battery.storage_change_j(),
                        inverter_storage_change_j,
                        -native.delta_h,
                    ]),
                    source_scale.max(inverter_scale),
                ),
                (sum([heat, -native.dissipated]), heat),
                (
                    sum([native.supplied, -shaft_work]),
                    sum(inverter.shaft_work_j.map(f64::abs)),
                ),
                (
                    step.balance_residual_j(),
                    source_scale.max(inverter_scale).max(heat),
                ),
            ];
            if [lead_current_a, lead_heat_j, inverter_storage_change_j]
                .iter()
                .chain(&step.discrete_bus_voltage_v)
                .chain(&step.heat_joules())
                .any(|v| !v.is_finite())
                || balances.iter().any(|&(residual, scale)| {
                    let allowance = self.controls.absolute_energy_tolerance_j
                        + self.controls.relative_energy_tolerance * scale;
                    !scale.is_finite()
                        || !allowance.is_finite()
                        || !residual.is_finite()
                        || residual.abs() > allowance
                })
            {
                return Err(invalid("battery, lead, inverter and shaft work do not close").into());
            }
            poll(&mut cancelled)?;
            Ok(step)
        })();
        match result {
            Ok(step) => {
                self.model.held = stamp;
                self.model.inverter.held = stamp.bridge;
                *next = candidate;
                Ok(step)
            }
            Err(error) => {
                let previous = self.model.held;
                self.model.apply_stamp(&previous);
                Err(error)
            }
        }
    }
}
