//! `motor_drive`: fs-system's [`CommutatedBatteryMotorDrive`], the
//! battery-fed commutated motor drive of fs-phs (an equivalent-circuit
//! battery and its lead, a two-leg converter of four switches and a DC motor
//! whose rotor and stator both turn) as one element with its two shafts and
//! its ten heat ports.
//!
//! The scene authors the drive, one gate program, the start and the owner's
//! controls; a step is the owner's whole program, so the scene's `step_s` is
//! its period. With `output_capacitance_f` the switches store charge (the
//! owner's eleven-state lane); without it they do not (nine states). Not
//! exposed: the gate lane, charge-curve storage, auxiliary loads, scheduled
//! programs, and resistance or switch laws: every resistance and switch
//! keeps its authored value.

use fs_ports::component::{Component, PortSpec};
use fs_system::battery::EquivalentCircuitBattery;
use fs_system::commutated_drive::{
    BodyDiode, CapacitiveCommutationParameters, CommutatedBatteryMotorDrive, CommutationInterval,
    CommutationParameters, DcMotor, DcMotorParameters, DeviceGateState, DeviceOperatingBounds,
    DeviceStorage, DriveParameters, GateProgram, HEAT_PORTS, NoDeviceStorage, PowerLoadControls,
    SwitchCycleControls, SwitchDevice,
};
use fs_system::{AdmissionError, NodeId, System, SystemBuilder};

use super::battery::{cell, halvings, solve, start as cell_start};
use super::{J, Refusal, fields, index_below, number, positive};

/// fs-system's test values for an interval's energy allowance when the
/// scene gives no `interval_tolerance`: absolute \[J\], relative.
const DEFAULT_INTERVAL: [f64; 2] = [1e-10, 1e-12];
/// fs-system's test values for the cycle when the scene gives no `cycle`.
const DEFAULT_CYCLE: SwitchCycleControls = SwitchCycleControls {
    max_intervals: 64,
    absolute_energy_tolerance_j: 1e-9,
    relative_energy_tolerance: 1e-11,
};
/// The owner's state without switch storage, in its order.
const NINE: [&str; 9] = [
    "remaining_charge_c",
    "polarization_charge_c",
    "battery_bus_charge_c",
    "converter_flux_wb",
    "input_bus_charge_c",
    "output_bus_charge_c",
    "motor_flux_wb",
    "rotor_momentum_n_m_s",
    "stator_momentum_n_m_s",
];
/// The switching-node charges the eleven-state lane adds after the bus
/// charges.
const NODES: [&str; 2] = ["input_node_charge_c", "output_node_charge_c"];
/// A drive's fields.
const FIELDS: [&str; 13] = [
    "kind",
    "battery",
    "lead_resistance_ohm",
    "converter",
    "motor",
    "output_capacitance_f",
    "program",
    "start",
    "energy_tolerance_j",
    "solve",
    "backtracking_halvings",
    "interval_tolerance",
    "cycle",
];
/// A drive's battery fields: the battery kind's curve and four constants.
const CELL: [&str; 6] = [
    "capacity_c",
    "soc_voltage",
    "series_resistance_ohm",
    "polarization_resistance_ohm",
    "polarization_capacitance_f",
    "bus_capacitance_f",
];

/// What a report and the scene's step admission read: the lane, the period,
/// the values that turn the owner's state into a state of charge, voltages,
/// currents and speeds, and the controls the run used.
#[derive(Clone, Debug)]
pub(super) struct Declared {
    eleven: bool,
    period_s: f64,
    capacity_c: f64,
    polarization_capacitance_f: f64,
    battery_bus_capacitance_f: f64,
    converter: CommutationParameters,
    motor: DcMotorParameters,
    electrical: PowerLoadControls,
    cycle: SwitchCycleControls,
    energy_tolerance_j: f64,
}

/// The state's names in the lane's order.
fn names(eleven: bool) -> Vec<&'static str> {
    let nodes: &[&str] = if eleven { &NODES } else { &[] };
    [&NINE[..6], nodes, &NINE[6..]].concat()
}

/// The object at `key`.
fn object<'j>(value: &'j J, key: &str, what: &str) -> Result<&'j J, String> {
    value
        .get(key)
        .filter(|object| object.as_object().is_some())
        .ok_or_else(|| format!("{what}.{key}: an object"))
}

/// `N` numbers at `key`.
fn numbers<const N: usize>(value: &J, key: &str, what: &str) -> Result<[f64; N], String> {
    value
        .get(key)
        .and_then(J::as_array)
        .and_then(|list| list.iter().map(J::as_f64).collect::<Option<Vec<_>>>())
        .and_then(|list| <[f64; N]>::try_from(list).ok())
        .ok_or_else(|| format!("{what}.{key}: {N} numbers"))
}

/// The owner's controls: the Newton budget and halvings of each interval's
/// solve with its energy allowance, and the cycle's; each object the scene
/// omits takes fs-system's test values, and one it gives has all its fields.
fn controls(
    definition: &J,
    what: &str,
) -> Result<(PowerLoadControls, SwitchCycleControls), String> {
    let [absolute, relative] = match definition.get("interval_tolerance") {
        None => DEFAULT_INTERVAL,
        Some(given) => {
            let what = format!("{what}.interval_tolerance");
            fields(given, &["absolute_j", "relative"], &what)?;
            [
                number(given, "absolute_j", &what)?,
                number(given, "relative", &what)?,
            ]
        }
    };
    let electrical = PowerLoadControls {
        solve: solve(definition, what)?,
        backtracking_halvings: halvings(definition, what)?,
        absolute_energy_tolerance_j: absolute,
        relative_energy_tolerance: relative,
    };
    let cycle = match definition.get("cycle") {
        None => DEFAULT_CYCLE,
        Some(given) => {
            let what = format!("{what}.cycle");
            fields(given, &["max_intervals", "absolute_j", "relative"], &what)?;
            SwitchCycleControls {
                max_intervals: given
                    .get("max_intervals")
                    .and_then(|count| index_below(count, 4097))
                    .filter(|&count| count > 0)
                    .ok_or_else(|| format!("{what}.max_intervals: a whole number, 1 to 4096"))?,
                absolute_energy_tolerance_j: number(given, "absolute_j", &what)?,
                relative_energy_tolerance: number(given, "relative", &what)?,
            }
        }
    };
    Ok((electrical, cycle))
}

/// One switch: its channel, its body diode and its limits.
fn switch(value: &J, what: &str) -> Result<SwitchDevice, String> {
    fields(
        value,
        &[
            "on_resistance_ohm",
            "off_resistance_ohm",
            "diode_forward_v",
            "diode_resistance_ohm",
            "max_reverse_v",
            "max_channel_a",
            "max_diode_a",
        ],
        what,
    )?;
    Ok(SwitchDevice {
        on_resistance_ohm: number(value, "on_resistance_ohm", what)?,
        off_resistance_ohm: number(value, "off_resistance_ohm", what)?,
        body_diode: BodyDiode {
            forward_voltage_v: number(value, "diode_forward_v", what)?,
            differential_resistance_ohm: number(value, "diode_resistance_ohm", what)?,
        },
        bounds: DeviceOperatingBounds {
            maximum_reverse_voltage_v: number(value, "max_reverse_v", what)?,
            maximum_channel_current_a: number(value, "max_channel_a", what)?,
            maximum_diode_current_a: number(value, "max_diode_a", what)?,
        },
    })
}

/// The converter: one switch for all four, or four in the owner's order
/// (input high, input low, output high, output low).
fn converter(definition: &J, what: &str) -> Result<CommutationParameters, String> {
    let value = object(definition, "converter", what)?;
    let what = format!("{what}.converter");
    fields(
        value,
        &[
            "inductance_h",
            "bus_capacitance_f",
            "winding_resistance_ohm",
            "switch",
        ],
        &what,
    )?;
    let given = value
        .get("switch")
        .ok_or_else(|| format!("{what}: missing switch"))?;
    let each = |index: usize, value: &J| switch(value, &format!("{what}.switch[{index}]"));
    let devices = match given.as_array() {
        None => [switch(given, &format!("{what}.switch"))?; 4],
        Some([a, b, c, d]) => [each(0, a)?, each(1, b)?, each(2, c)?, each(3, d)?],
        Some(_) => return Err(format!("{what}.switch: one switch, or a list of four")),
    };
    Ok(CommutationParameters {
        inductance_h: number(value, "inductance_h", &what)?,
        bus_capacitance_f: numbers(value, "bus_capacitance_f", &what)?,
        winding_resistance_ohm: number(value, "winding_resistance_ohm", &what)?,
        devices,
    })
}

/// The motor's constants.
fn motor(definition: &J, what: &str) -> Result<DcMotorParameters, String> {
    let value = object(definition, "motor", what)?;
    let what = format!("{what}.motor");
    let keys = [
        "inductance_h",
        "winding_resistance_ohm",
        "rotor_inertia_kg_m2",
        "stator_inertia_kg_m2",
        "coupling_nm_per_a",
        "bearing_drag_nm_s_per_rad",
    ];
    fields(value, &keys, &what)?;
    Ok(DcMotorParameters {
        inductance_h: number(value, keys[0], &what)?,
        winding_resistance_ohm: number(value, keys[1], &what)?,
        rotor_inertia_kg_m2: number(value, keys[2], &what)?,
        stator_inertia_kg_m2: number(value, keys[3], &what)?,
        coupling_nm_per_a: number(value, keys[4], &what)?,
        bearing_drag_nm_s_per_rad: number(value, keys[5], &what)?,
    })
}

/// The gate program, every step: its intervals, each with the four gate
/// commands in the owner's switch order, and the owner's substeps.
fn program(definition: &J, what: &str) -> Result<GateProgram, String> {
    let value = object(definition, "program", what)?;
    let what = format!("{what}.program");
    fields(value, &["intervals", "substeps"], &what)?;
    let rows = value
        .get("intervals")
        .and_then(J::as_array)
        .filter(|rows| !rows.is_empty())
        .ok_or_else(|| format!("{what}.intervals: a nonempty list"))?;
    let intervals = rows
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let what = format!("{what}.intervals[{index}]");
            fields(row, &["on", "duration_s"], &what)?;
            let on = row
                .get("on")
                .and_then(J::as_array)
                .and_then(|flags| {
                    let flags = flags.iter().map(|flag| match flag {
                        J::Bool(on) => Some(*on),
                        _ => None,
                    });
                    <[bool; 4]>::try_from(flags.collect::<Option<Vec<_>>>()?).ok()
                })
                .ok_or_else(|| {
                    format!(
                        "{what}.on: four booleans, input high, input low, output high, output low"
                    )
                })?;
            Ok(CommutationInterval {
                gates: DeviceGateState { on },
                duration_s: positive(row, "duration_s", &what)?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let substeps_per_interval = value
        .get("substeps")
        .and_then(|substeps| index_below(substeps, 1025))
        .filter(|&substeps| substeps > 0)
        .ok_or_else(|| format!("{what}.substeps: a whole number, 1 to 1024"))?;
    Ok(GateProgram {
        auxiliary_loads: vec![None; intervals.len()],
        intervals,
        substeps_per_interval,
    })
}

/// The start at `value`, in the owner's order: the battery's three charges,
/// the converter's flux and bus charges (and with `storage` its switching
/// nodes' charges), the motor's flux and its two momenta, each mapped by its
/// owner's `state` from the authored state of charge, voltages, currents and
/// shaft speeds.
fn initial(
    value: &J,
    battery: &EquivalentCircuitBattery,
    converter: CommutationParameters,
    storage: Option<[f64; 4]>,
    motor: DcMotorParameters,
    what: &str,
) -> Result<Vec<f64>, String> {
    let mut allowed = vec![
        "soc",
        "polarization_v",
        "bus_v",
        "converter_current_a",
        "converter_bus_v",
        "motor_current_a",
        "shaft_rad_s",
    ];
    if storage.is_some() {
        allowed.push("node_v");
    }
    fields(value, &allowed, what)?;
    let refused = |e| format!("{what}: {e}");
    let mut state = cell_start(battery, value, what)?.to_vec();
    let current_a = number(value, "converter_current_a", what)?;
    let bus_v = numbers(value, "converter_bus_v", what)?;
    if let Some(output_capacitance_f) = storage {
        let switches = CapacitiveCommutationParameters {
            converter,
            output_capacitance_f,
        };
        let node_v = numbers(value, "node_v", what)?;
        state.extend(switches.state(current_a, bus_v, node_v).map_err(refused)?);
    } else {
        state.extend(converter.state(current_a, bus_v).map_err(refused)?);
    }
    let motor_current_a = number(value, "motor_current_a", what)?;
    let shaft_rad_s = numbers(value, "shaft_rad_s", what)?;
    let motor = DcMotor::new(motor).map_err(refused)?;
    state.extend(motor.state(motor_current_a, shaft_rad_s).map_err(refused)?);
    Ok(state)
}

/// A `motor_drive` added to `builder` as `name` at its start, with its ports
/// and what its reports read.
pub(super) fn add(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>, Declared), Refusal> {
    fields(definition, &FIELDS, what)?;
    let energy_tolerance_j = positive(definition, "energy_tolerance_j", what)?;
    let (electrical, cycle) = controls(definition, what)?;
    let value = object(definition, "battery", what)?;
    let battery_what = format!("{what}.battery");
    fields(value, &CELL, &battery_what)?;
    let (battery, cell) = cell(
        value,
        &battery_what,
        electrical.solve,
        electrical.backtracking_halvings,
        energy_tolerance_j,
    )?;
    let converter = converter(definition, what)?;
    let motor = motor(definition, what)?;
    let storage = match definition.get("output_capacitance_f") {
        None => None,
        Some(_) => Some(numbers::<4>(definition, "output_capacitance_f", what)?),
    };
    let program = program(definition, what)?;
    let start = object(definition, "start", what)?;
    let state = initial(
        start,
        &battery,
        converter,
        storage,
        motor,
        &format!("{what}.start"),
    )?;
    let mut declared = Declared {
        eleven: storage.is_some(),
        period_s: 0.0,
        capacity_c: cell.capacity_c,
        polarization_capacitance_f: cell.parameters.polarization_capacitance_f,
        battery_bus_capacitance_f: cell.parameters.bus_capacitance_f,
        converter,
        motor,
        electrical,
        cycle,
        energy_tolerance_j,
    };
    let parameters = DriveParameters {
        capacity_c: cell.capacity_c,
        soc_voltage: cell.soc_voltage,
        battery: cell.parameters,
        lead_resistance_ohm: number(definition, "lead_resistance_ohm", what)?,
        converter,
        motor,
    };
    let refused = |e| format!("{what}: {e}");
    let admission = |e: AdmissionError| Refusal::Admission(e.to_string());
    let shape = |_| format!("{what}: the start does not fit the lane");
    let (node, ports, period_s) = if let Some(output_capacitance_f) = storage {
        let storage = DeviceStorage::Capacitance(output_capacitance_f);
        let drive = CommutatedBatteryMotorDrive::new(
            parameters,
            storage,
            program,
            electrical,
            cycle,
            energy_tolerance_j,
        )
        .map_err(refused)?;
        let (ports, period_s) = (drive.ports().to_vec(), drive.period_s());
        let state = <[f64; 11]>::try_from(state).map_err(shape)?;
        let node = builder.add(name, drive, state).map_err(admission)?;
        (node, ports, period_s)
    } else {
        let drive = CommutatedBatteryMotorDrive::new(
            parameters,
            NoDeviceStorage,
            program,
            electrical,
            cycle,
            energy_tolerance_j,
        )
        .map_err(refused)?;
        let (ports, period_s) = (drive.ports().to_vec(), drive.period_s());
        let state = <[f64; 9]>::try_from(state).map_err(shape)?;
        let node = builder.add(name, drive, state).map_err(admission)?;
        (node, ports, period_s)
    };
    declared.period_s = period_s;
    Ok((node, ports, declared))
}

/// Refuses a scene whose `step_s` is not the drive's period within the
/// owner's 64 rounding units, or whose `duration_s` is not a whole number of
/// steps within 1e-13 of itself (inside the run's 1e-12 end slack): a step of
/// the drive is its whole program, which the owner never cuts.
pub(super) fn admit(
    declared: &Declared,
    step_s: f64,
    duration_s: f64,
    what: &str,
) -> Result<(), String> {
    let period_s = declared.period_s;
    let steps = (duration_s / step_s).round();
    if (step_s - period_s).abs() > 64.0 * f64::EPSILON * period_s
        || (duration_s - steps * step_s).abs() > 1e-13 * duration_s
    {
        return Err(format!(
            "{what}: a step is one period of its program, {period_s} s: step_s must be the \
             period and duration_s a whole number of steps"
        ));
    }
    Ok(())
}

/// The state by name, in the owner's order; the state of charge, the
/// battery's polarization and bus voltages, the converter's and the motor's
/// currents, both shaft speeds and their difference, and without switch
/// storage the two bus voltages; each heat port's last heat rate (null
/// before the first step); the controls the run used and the period.
pub(super) fn values(system: &System, node: NodeId, declared: &Declared) -> Vec<(String, f64)> {
    let names = names(declared.eleven);
    let state = if declared.eleven {
        system.state::<[f64; 11]>(node).map(|state| state.to_vec())
    } else {
        system.state::<[f64; 9]>(node).map(|state| state.to_vec())
    }
    .unwrap_or_else(|| vec![f64::NAN; names.len()]);
    let at = |name: &str| {
        names
            .iter()
            .position(|known| *known == name)
            .map_or(f64::NAN, |index| state[index])
    };
    let (motor, converter) = (declared.motor, declared.converter);
    let rotor = at("rotor_momentum_n_m_s") / motor.rotor_inertia_kg_m2;
    let stator = at("stator_momentum_n_m_s") / motor.stator_inertia_kg_m2;
    let mut values: Vec<(String, f64)> = names
        .iter()
        .zip(&state)
        .map(|(name, value)| ((*name).to_owned(), *value))
        .collect();
    let mut derived = vec![
        ("soc", at("remaining_charge_c") / declared.capacity_c),
        (
            "polarization_v",
            at("polarization_charge_c") / declared.polarization_capacitance_f,
        ),
        (
            "battery_bus_v",
            at("battery_bus_charge_c") / declared.battery_bus_capacitance_f,
        ),
        (
            "converter_current_a",
            at("converter_flux_wb") / converter.inductance_h,
        ),
        ("motor_current_a", at("motor_flux_wb") / motor.inductance_h),
        ("rotor_rad_s", rotor),
        ("stator_rad_s", stator),
        ("relative_rad_s", rotor - stator),
    ];
    if !declared.eleven {
        // Without switch storage each bus charge is its capacitor's alone.
        let [input, output] = converter.bus_capacitance_f;
        derived.push(("input_bus_v", at("input_bus_charge_c") / input));
        derived.push(("output_bus_v", at("output_bus_charge_c") / output));
    }
    values.extend(
        derived
            .into_iter()
            .map(|(name, value)| (name.to_owned(), value)),
    );
    for port in HEAT_PORTS {
        let rate = system
            .port(node, port)
            .filter(|_| system.steps() > 0)
            .and_then(|port| system.last_output(port).first().copied())
            .map_or(f64::NAN, |output| -output);
        values.push((format!("{port}_w"), rate));
    }
    // The electrical controls, their Newton budget and the cycle's.
    let (e, c) = (declared.electrical, declared.cycle);
    let s = e.solve;
    let controls = [
        ("solve_relative", s.relative),
        ("solve_absolute", s.absolute),
        ("solve_max_iterations", s.max_iterations as f64),
        ("backtracking_halvings", f64::from(e.backtracking_halvings)),
        ("interval_absolute_j", e.absolute_energy_tolerance_j),
        ("interval_relative", e.relative_energy_tolerance),
        ("cycle_max_intervals", c.max_intervals as f64),
        ("cycle_absolute_j", c.absolute_energy_tolerance_j),
        ("cycle_relative", c.relative_energy_tolerance),
        ("energy_tolerance_j", declared.energy_tolerance_j),
        ("period_s", declared.period_s),
    ];
    values.extend(controls.map(|(name, value)| (name.to_owned(), value)));
    values
}
