//! `battery`: fs-system's [`EquivalentCircuitBattery`], the OCV/SOC curve of
//! `fs_phs::circuit::battery` with a series resistance, one polarization RC
//! pair and a bus capacitor, as one element with its terminal and its two
//! heat ports.
//!
//! The scene authors the curve, the four constants, the start and the
//! owner's controls; every step is the owner's own. The resistances are
//! constant: the Component's temperature laws are not exposed. The controls
//! a scene omits take fs-system's own test values, and every report names
//! the controls the run used.

use fs_system::battery::{BatteryParameters, EquivalentCircuitBattery, SolveBudget};
use fs_system::{NodeId, System};

use super::{J, Refusal, fields, number, number_or, positive};

/// The owner's Newton budget when the scene gives no `solve`.
const DEFAULT_SOLVE: SolveBudget = SolveBudget {
    relative: 1e-13,
    absolute: 1e-12,
    max_iterations: 40,
};
/// The owner's step halvings when the scene gives no `backtracking_halvings`.
const DEFAULT_HALVINGS: u8 = 10;
/// The owner's state, in its order (the remaining chemical charge, the
/// polarization charge and the bus charge \[C\]), as reported and as the
/// battery's adaptive accounts.
pub(super) const CHARGES: [&str; 3] = [
    "remaining_charge_c",
    "polarization_charge_c",
    "bus_charge_c",
];

/// What a report reads besides the state: the capacity and the two
/// capacitances that turn the owner's charges into a state of charge and
/// voltages, and the controls the run used.
#[derive(Clone, Copy, Debug)]
pub(super) struct Declared {
    capacity_c: f64,
    polarization_capacitance_f: f64,
    bus_capacitance_f: f64,
    solve: SolveBudget,
    halvings: u8,
    energy_tolerance_j: f64,
}

/// A cell's capacity, curve and constants, which a `motor_drive` also reads.
pub(super) struct Cell {
    pub(super) capacity_c: f64,
    pub(super) soc_voltage: Vec<(f64, f64)>,
    pub(super) parameters: BatteryParameters,
}

/// A `battery` and its start, `[remaining charge, polarization charge, bus
/// charge]` \[C\], from its state of charge and two voltages.
pub(super) fn battery(
    definition: &J,
    what: &str,
) -> Result<(EquivalentCircuitBattery, [f64; 3], Declared), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "capacity_c",
            "soc_voltage",
            "series_resistance_ohm",
            "polarization_resistance_ohm",
            "polarization_capacitance_f",
            "bus_capacitance_f",
            "soc",
            "polarization_v",
            "bus_v",
            "energy_tolerance_j",
            "solve",
            "backtracking_halvings",
        ],
        what,
    )?;
    let solve = solve(definition, what)?;
    let halvings = halvings(definition, what)?;
    let energy_tolerance_j = positive(definition, "energy_tolerance_j", what)?;
    let (battery, cell) = cell(definition, what, solve, halvings, energy_tolerance_j)?;
    let state = start(&battery, definition, what)?;
    let declared = Declared {
        capacity_c: cell.capacity_c,
        polarization_capacitance_f: cell.parameters.polarization_capacitance_f,
        bus_capacitance_f: cell.parameters.bus_capacitance_f,
        solve,
        halvings,
        energy_tolerance_j,
    };
    Ok((battery, state, declared))
}

/// The cell at `definition` (its curve and four constants, whose fields the
/// caller has checked), built under the owner's controls.
pub(super) fn cell(
    definition: &J,
    what: &str,
    solve: SolveBudget,
    halvings: u8,
    energy_tolerance_j: f64,
) -> Result<(EquivalentCircuitBattery, Cell), String> {
    let capacity_c = positive(definition, "capacity_c", what)?;
    let soc_voltage = definition
        .get("soc_voltage")
        .and_then(J::as_array)
        .and_then(|knots| {
            knots
                .iter()
                .map(|knot| match knot.as_array()? {
                    [soc, volts] => Some((soc.as_f64()?, volts.as_f64()?)),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
        })
        .ok_or_else(|| {
            format!("{what}.soc_voltage: [[soc, V], ...] pairs of numbers, soc from 0 to 1")
        })?;
    let parameters = BatteryParameters {
        series_resistance_ohm: number(definition, "series_resistance_ohm", what)?,
        polarization_resistance_ohm: number(definition, "polarization_resistance_ohm", what)?,
        polarization_capacitance_f: number(definition, "polarization_capacitance_f", what)?,
        bus_capacitance_f: number(definition, "bus_capacitance_f", what)?,
    };
    let battery = EquivalentCircuitBattery::new(
        capacity_c,
        &soc_voltage,
        parameters,
        solve,
        halvings,
        energy_tolerance_j,
    )
    .map_err(|e| format!("{what}: {e}"))?;
    let cell = Cell {
        capacity_c,
        soc_voltage,
        parameters,
    };
    Ok((battery, cell))
}

/// The owner's Newton budget at `solve`, or fs-system's test value.
pub(super) fn solve(definition: &J, what: &str) -> Result<SolveBudget, String> {
    let solve = match definition.get("solve") {
        None => DEFAULT_SOLVE,
        Some(budget) => {
            let what = format!("{what}.solve");
            fields(budget, &["relative", "absolute", "max_iterations"], &what)?;
            SolveBudget {
                relative: number(budget, "relative", &what)?,
                absolute: number(budget, "absolute", &what)?,
                max_iterations: budget
                    .get("max_iterations")
                    .and_then(J::as_f64)
                    .filter(|n| n.fract() == 0.0 && (1.0..=1000.0).contains(n))
                    .map(|n| n as usize)
                    .ok_or_else(|| format!("{what}.max_iterations: a whole number, 1 to 1000"))?,
            }
        }
    };
    Ok(solve)
}

/// The owner's step halvings at `backtracking_halvings`, or fs-system's test
/// value.
pub(super) fn halvings(definition: &J, what: &str) -> Result<u8, String> {
    let halvings = match definition.get("backtracking_halvings") {
        None => DEFAULT_HALVINGS,
        Some(value) => value
            .as_f64()
            .filter(|n| n.fract() == 0.0 && (0.0..=64.0).contains(n))
            .map(|n| n as u8)
            .ok_or_else(|| format!("{what}.backtracking_halvings: a whole number, 0 to 64"))?,
    };
    Ok(halvings)
}

/// The start at `definition`, the three charges mapped by the owner's
/// `state` from the state of charge and two voltages.
pub(super) fn start(
    battery: &EquivalentCircuitBattery,
    definition: &J,
    what: &str,
) -> Result<[f64; 3], String> {
    let state = battery
        .state(
            number(definition, "soc", what)?,
            number_or(definition, "polarization_v", 0.0, what)?,
            number(definition, "bus_v", what)?,
        )
        .map_err(|e| format!("{what}: {e}"))?;
    Ok(state)
}

/// The three charges, in the owner's order; then the state of charge and the
/// two capacitor voltages; the last step's terminal current (positive charges
/// the battery), its mean bus voltage and the two heat rates (null before the
/// first step); then the controls the run used.
pub(super) fn values(system: &System, node: NodeId, declared: &Declared) -> Vec<(String, f64)> {
    let state = system
        .state::<[f64; 3]>(node)
        .copied()
        .unwrap_or([f64::NAN; 3]);
    let last = |port: &str, input: bool| {
        system
            .port(node, port)
            .filter(|_| system.steps() > 0)
            .and_then(|port| {
                let values = if input {
                    system.last_input(port)
                } else {
                    system.last_output(port)
                };
                values.first().copied()
            })
            .unwrap_or(f64::NAN)
    };
    let derived = [
        ("soc", state[0] / declared.capacity_c),
        (
            "polarization_v",
            state[1] / declared.polarization_capacitance_f,
        ),
        ("bus_v", state[2] / declared.bus_capacitance_f),
        ("current_a", last("terminal", true)),
        ("terminal_v", last("terminal", false)),
        ("heat_series_w", -last("heat_series", false)),
        ("heat_polarization_w", -last("heat_polarization", false)),
        ("solve_relative", declared.solve.relative),
        ("solve_absolute", declared.solve.absolute),
        ("solve_max_iterations", declared.solve.max_iterations as f64),
        ("backtracking_halvings", f64::from(declared.halvings)),
        ("energy_tolerance_j", declared.energy_tolerance_j),
    ];
    CHARGES
        .into_iter()
        .zip(state)
        .chain(derived)
        .map(|(name, value)| (name.to_owned(), value))
        .collect()
}
