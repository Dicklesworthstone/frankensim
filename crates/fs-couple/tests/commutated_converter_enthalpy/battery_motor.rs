//! Finite chemical/shaft power, diode heat and the original finite caloric owners.
#[path = "battery_motor/capacitance.rs"]
mod capacitance;
#[path = "battery_motor/charge.rs"]
mod charge;
#[path = "battery_motor/gate.rs"]
mod gate;

use super::*;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
            PreparedHeatNetwork,
        },
        material::ElectricalContactNetwork,
    },
    converter::commutation::battery_motor::*,
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_phs::circuit::{
    battery::{BatteryOcvCurve, BatteryParameters, BatterySource},
    converter::{
        AveragedBuckBoost, BridgeModulation, BuckBoostParameters,
        battery_fed::motor::BatteryMotorDrive, switching::battery_motor::BatteryMotorSwitchCommand,
    },
    motor::{DcMotor, DcMotorParameters},
};
use std::collections::BTreeMap;

fn motor() -> DcMotorParameters {
    DcMotorParameters {
        inductance_h: 0.05,
        winding_resistance_ohm: 0.1,
        rotor_inertia_kg_m2: 0.4,
        stator_inertia_kg_m2: 0.8,
        coupling_nm_per_a: 0.2,
        bearing_drag_nm_s_per_rad: 0.0003,
    }
}
fn drive_resistors(f: &Fixture) -> [ConverterThermalChannel<'_>; 5] {
    // Explicit shared bodies: battery series/lead, polarization/motor winding.
    [IDS[0], IDS[5], IDS[6], IDS[5], IDS[6]].map(|cell| ConverterThermalChannel {
        cell,
        resistance: ConverterResistance::Resistance {
            key: &f.keys[0],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    })
}
fn drive_controls(cells: usize) -> CircuitEnthalpySolve {
    let mut c = controls();
    c.maximum_cells = cells;
    c.newton.linear_restart = cells;
    c.maximum_circuit_trials = 10_000;
    c
}

struct SteadyDrive {
    native: BatteryMotorDrive,
    state: [f64; 9],
    command: BatteryMotorSwitchCommand,
    gates: DeviceGateState,
    cell_current_a: f64,
    heat_w: [f64; 14],
}

// G1: solve original KCL/KVL algebra independently, with explicit OFF leakage.
// Motoring uses an output both-OFF interval and its high body diode; regeneration
// exchanges the roles. Both external shaft reactions are authored, retaining
// finite absolute rotor/stator speeds. Only chemical charge changes in time.
fn steady_drive(regeneration: bool) -> SteadyDrive {
    let current = if regeneration { -1.0 } else { 1.0 };
    let (g_off, g_diode, r_on, winding_r) = (0.01, 5.0, 0.1, 0.1);
    let g = g_off + g_diode;
    let source_r = 0.3; // Battery series, polarization and physical lead.
    let left = if regeneration {
        (g * 12.0 + g_diode * 0.5 - current * (1.0 + g * source_r))
            / (g + g_off * (1.0 + g * source_r))
    } else {
        (12.0 - (source_r + r_on) * current) / (1.0 + g_off * (source_r + r_on))
    };
    let cell_current = current + g_off * left;
    let input = 12.0 - source_r * cell_current;
    let battery_bus = 12.0 - 0.2 * cell_current;
    let right = left - winding_r * current;
    let output = if regeneration {
        right + r_on * (-current + g_off * right)
    } else {
        ((2.0 * g_off + g_diode) * right - g_diode * 0.5 - current) / g
    };
    let auxiliary_current = 0.1;
    let motor_current = current - g_off * right - auxiliary_current;
    let m = motor();
    let relative = (output - m.winding_resistance_ohm * motor_current) / m.coupling_nm_per_a;
    let torque = m.bearing_drag_nm_s_per_rad * relative - m.coupling_nm_per_a * motor_current;
    let gates = DeviceGateState {
        on: if regeneration {
            [false, false, true, false]
        } else {
            [true, false, false, false]
        },
    };
    let p = parameters();
    let battery = BatterySource::new(
        BatteryOcvCurve::new(100.0, &[(0.0, 12.0), (1.0, 12.0)]).unwrap(),
        BatteryParameters {
            series_resistance_ohm: 0.1,
            polarization_resistance_ohm: 0.1,
            polarization_capacitance_f: 2.0,
            bus_capacitance_f: 0.5,
        },
    )
    .unwrap();
    let electrical = AveragedBuckBoost::new(
        BuckBoostParameters {
            inductance_h: p.inductance_h,
            bus_capacitance_f: p.bus_capacitance_f,
            winding_resistance_ohm: p.winding_resistance_ohm,
            input_on_resistance_ohm: [0.1; 2],
            output_on_resistance_ohm: [0.1; 2],
        },
        BridgeModulation {
            input_high_fraction: 0.0,
            output_high_fraction: 0.0,
        },
    )
    .unwrap()
    .with_battery_source(battery, 0.1)
    .unwrap();
    let electrical_state = electrical
        .state(
            0.8,
            0.1 * cell_current,
            battery_bus,
            current,
            [input, output],
        )
        .unwrap();
    let machine = DcMotor::new(m).unwrap();
    let motor_state = machine.state(motor_current, [relative + 3.0, 3.0]).unwrap();
    let native = electrical.with_dc_motor(machine).unwrap();
    let state = native.state(&electrical_state, &motor_state).unwrap();
    let vds = [input - left, left, output - right, right];
    let mut heat = [0.0; 14];
    heat[0] = winding_r * current * current;
    for i in 0..4 {
        heat[1 + i] = vds[i] * vds[i] / if gates.on[i] { 0.1 } else { 100.0 };
        heat[10 + i] = -vds[i] * ((-vds[i] - 0.5) / 0.2).max(0.0);
    }
    heat[5] = 0.1 * cell_current * cell_current;
    heat[6] = heat[5];
    heat[7] = heat[5];
    heat[8] = 0.1 * motor_current * motor_current;
    heat[9] = m.bearing_drag_nm_s_per_rad * relative * relative;
    SteadyDrive {
        native,
        state,
        gates,
        cell_current_a: cell_current,
        heat_w: heat,
        command: BatteryMotorSwitchCommand {
            shaft_torque_nm: [torque, -torque],
            auxiliary_load: Some(
                RegulatedPowerLoad::try_new(output * auxiliary_current, 10.0, 0.25, 100.0).unwrap(),
            ),
        },
    }
}
fn per_body(heat: [f64; 14], routed_load: f64) -> [f64; 7] {
    [
        heat[0] + heat[9],
        heat[1] + heat[10],
        heat[2] + heat[11],
        heat[3] + heat[12] + routed_load,
        heat[4] + heat[13],
        heat[5] + heat[7],
        heat[6] + heat[8],
    ]
}

#[test]
fn independent_motoring_and_regenerative_dead_time_calorimetry_closes_all_fourteen_heats() {
    for regeneration in [false, true] {
        let f = fixture(1.0, 0.0, 350.0, 100.0);
        let thermal = f.thermal();
        let context = common::point(300.0);
        let cells = f.cells(&context);
        let resistors = drive_resistors(&f);
        let devices = f.devices();
        let network = CommutatedBatteryMotorThermalNetwork {
            cells: &cells,
            links: &[],
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let steady = steady_drive(regeneration);
        let drive = network
            .bind(&thermal, &steady.native, steady.state, || false)
            .unwrap();
        let dt = 1.0 / 64.0;
        let intervals = [CommutationInterval {
            gates: steady.gates,
            duration_s: dt,
        }];
        let run = |receiver| {
            network
                .step_candidate(
                    &thermal,
                    dt,
                    &drive,
                    &intervals,
                    &[steady.command],
                    4,
                    receiver,
                    power_controls(),
                    cycle_controls(),
                    drive_controls(7),
                    &[0.0; 7],
                    || false,
                )
                .unwrap()
        };
        let exported = run(None);
        let heated = run(Some(IDS[3]));
        let step = &heated.report().electrical;
        let heat = steady.heat_w.map(|power| power * dt);
        for (actual, expected) in step.heat_joules().into_iter().zip(heat) {
            near(actual, expected);
        }
        let load = steady.command.auxiliary_load.unwrap().requested_power_w() * dt;
        near(step.load_work_j, load);
        let by_body = per_body(heat, load);
        for i in 0..7 {
            near(
                f.cp * (heated.thermal().temperature_k(IDS[i]).unwrap() - 300.0),
                by_body[i],
            );
        }
        near(
            heated.drive().state()[0],
            steady.state[0] - dt * steady.cell_current_a,
        );
        near(step.storage_change_j[0], -12.0 * dt * steady.cell_current_a);
        for i in 1..9 {
            near(heated.drive().state()[i], steady.state[i]);
        }
        assert_eq!(step.storage_change_j[0] > 0.0, regeneration);
        assert_eq!(step.motor_electrical_work_j < 0.0, regeneration);
        assert!(if regeneration {
            heat[10] > 0.0
        } else {
            heat[12] > 0.0
        });
        let m = motor();
        near(
            step.shaft_work_j[0],
            dt * steady.command.shaft_torque_nm[0] * steady.state[7] / m.rotor_inertia_kg_m2,
        );
        near(
            step.shaft_work_j[1],
            dt * steady.command.shaft_torque_nm[1] * steady.state[8] / m.stator_inertia_kg_m2,
        );
        near(
            step.native.supplied,
            step.shaft_work_j.iter().sum::<f64>() - load,
        );
        near(
            heated.report().deposited_joule_j,
            heat.iter().sum::<f64>() + load,
        );
        near(
            heated.report().deposited_joule_j - exported.report().deposited_joule_j,
            load,
        );
        near(
            heated.thermal().energy_j() - thermal.energy_j(),
            heated.report().deposited_joule_j,
        );
        near(heated.report().total_energy_residual_j, 0.0);
        assert_eq!(heated.report().circuit_trials % 4, 0);
        assert_eq!(heated.drive().battery_curve().capacity_c(), 100.0);
    }
}

#[test]
fn diode_threshold_feedback_reaches_native_motion_when_every_resistance_stays_constant() {
    let f = fixture(0.001, 0.02, 350.0, 100.0);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f);
    let devices = f.devices();
    let network = CommutatedBatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
        bearing_cell: IDS[0],
    };
    let steady = steady_drive(false);
    let drive = network
        .bind(&thermal, &steady.native, steady.state, || false)
        .unwrap();
    let dt = 1.0 / 1024.0;
    let result = network
        .step_candidate(
            &thermal,
            dt,
            &drive,
            &[CommutationInterval {
                gates: steady.gates,
                duration_s: dt,
            }],
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(7),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    let step = &result.report().electrical;
    let next = result.drive().state();
    let midpoint = |i: usize| f64::midpoint(steady.state[i], next[i]);
    let p = result.drive().parameters();
    let b = result.drive().battery_parameters();
    let m = result.drive().motor_parameters();
    let current = midpoint(3) / p.inductance_h;
    let input = midpoint(4) / p.bus_capacitance_f[0];
    let output = midpoint(5) / p.bus_capacitance_f[1];
    let vf = p.devices[2].body_diode.forward_voltage_v;
    let left = (10.0 * input - current) / 10.01;
    let right = (5.01 * output + 5.0 * vf + current) / 5.02;
    let vds = [input - left, left, output - right, right];
    let mut expected = [0.0; 14];
    expected[0] = dt * 0.1 * current * current;
    for i in 0..4 {
        let t = result.thermal().temperature_k(IDS[i + 1]).unwrap();
        near(p.devices[i].on_resistance_ohm, 0.1);
        near(p.devices[i].off_resistance_ohm, 100.0);
        near(
            p.devices[i].body_diode.forward_voltage_v,
            0.5 + 0.02 * (t - 300.0),
        );
        let diode_current = ((-vds[i] - p.devices[i].body_diode.forward_voltage_v) / 0.2).max(0.0);
        expected[1 + i] = dt * vds[i] * vds[i] / if steady.gates.on[i] { 0.1 } else { 100.0 };
        expected[10 + i] = dt * -vds[i] * diode_current;
        near(step.last_devices[i].diode_current_a, diode_current);
        assert_eq!(
            drive.devices()[i].source().card_identity(),
            result.drive().devices()[i].source().card_identity()
        );
    }
    let pol = midpoint(1) / b.polarization_capacitance_f;
    let battery_bus = midpoint(2) / b.bus_capacitance_f;
    let cell_current = (12.0 - pol - battery_bus) / 0.1;
    let lead_current = (battery_bus - input) / 0.1;
    let motor_current = midpoint(6) / m.inductance_h;
    let relative = midpoint(7) / m.rotor_inertia_kg_m2 - midpoint(8) / m.stator_inertia_kg_m2;
    expected[5] = dt * 0.1 * cell_current * cell_current;
    expected[6] = dt * pol * pol / 0.1;
    expected[7] = dt * 0.1 * lead_current * lead_current;
    expected[8] = dt * 0.1 * motor_current * motor_current;
    expected[9] = dt * m.bearing_drag_nm_s_per_rad * relative * relative;
    for (actual, expected) in step.heat_joules().into_iter().zip(expected) {
        near(actual, expected);
    }
    for (i, q) in per_body(expected, 0.0).into_iter().enumerate() {
        near(
            f.cp * (result.thermal().temperature_k(IDS[i]).unwrap() - 300.0),
            q,
        );
    }
    for r in result.drive().resistors() {
        near(r.resistance_ohm(), 0.1);
    }
    assert!(result.thermal().temperature_k(IDS[3]).unwrap() > 300.1);
    assert!(
        (next[3] - steady.state[3]).abs() > 1e-8,
        "threshold feedback must change actual converter flux"
    );
    near(result.report().total_energy_residual_j, 0.0);
}

#[test]
fn device_domain_limits_clock_budget_and_final_cancellation_publish_no_history() {
    for (upper, limit) in [(300.0001, 100.0), (350.0, 0.01), (350.0, 100.0)] {
        let f = fixture(0.001, 0.02, upper, limit);
        let thermal = f.thermal();
        let before = thermal.cells().clone();
        let context = common::point(300.0);
        let cells = f.cells(&context);
        let resistors = drive_resistors(&f);
        let devices = f.devices();
        let network = CommutatedBatteryMotorThermalNetwork {
            cells: &cells,
            links: &[],
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let steady = steady_drive(false);
        let drive = network
            .bind(&thermal, &steady.native, steady.state, || false)
            .unwrap();
        let dt = 1.0 / 1024.0;
        let intervals = [CommutationInterval {
            gates: steady.gates,
            duration_s: dt,
        }];
        let run = |target, controls, cancel: &mut dyn FnMut() -> bool| {
            network.step_candidate(
                &thermal,
                target,
                &drive,
                &intervals,
                &[steady.command],
                4,
                None,
                power_controls(),
                cycle_controls(),
                controls,
                &[0.0; 7],
                cancel,
            )
        };
        if upper < 301.0 || limit < 1.0 {
            assert!(run(dt, drive_controls(7), &mut || false).is_err());
        } else {
            let mut polls = 0;
            let accepted = run(dt, drive_controls(7), &mut || {
                polls += 1;
                false
            })
            .unwrap();
            assert!(run(dt.next_up(), drive_controls(7), &mut || false).is_err());
            let mut short = drive_controls(7);
            short.maximum_circuit_trials = 3;
            assert!(matches!(
                run(dt, short, &mut || false),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            let mut count = 0;
            assert!(matches!(
                run(dt, drive_controls(7), &mut || {
                    count += 1;
                    count >= polls
                }),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let retry = run(dt, drive_controls(7), &mut || false).unwrap();
            assert_eq!(retry.drive().state(), accepted.drive().state());
            assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
        }
        assert_eq!(drive.state(), &steady.state);
        assert_eq!(drive.time_s(), 0.0);
        assert_eq!(thermal.cells(), &before);
    }
}

// A fixed, authored finite thermal link through the ORIGINAL HeatAssembly seam.
// The caloric owner resolves both inventories; this provider integrates nothing.
struct FixedCooling<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    conductance_w_k: f64,
}
impl HeatAssembly for FixedCooling<'_, '_, '_> {
    type Sample = HeatNetworkSample;
    type Error = HeatNetworkError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.prepared
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<HeatNetworkSample>, CircuitEnthalpyError<HeatNetworkError>>
    {
        self.prepared
            .resolve_with_link_conductances(
                cx,
                thermal,
                &BTreeMap::from([(77, self.conductance_w_k)]),
            )
            .map_err(|error| match error {
                HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
                other => CircuitEnthalpyError::Source(other),
            })
    }
}

#[test]
fn commutated_battery_drive_and_finite_liquid_coolant_share_the_original_heat_assembly() {
    const COOLANT: CellId = CellId(80);
    let f = fixture(1.0, 0.0, 350.0, 100.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in f.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let fluid = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-coolant".into(),
            phase: "liquid".into(),
            process: "authored finite heat-capacity fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let fluid_curve = ResolvedHeatCapacityCurve::try_new(
        &fluid,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-6,
    )
    .unwrap();
    let context = common::point(300.0);
    let fluid_context = common::point(299.5);
    let mut cells = f.cells(&context).into_iter().collect::<Vec<_>>();
    cells.push(MaterialCell {
        cell: COOLANT,
        material_id: MaterialId(1),
        curve: fluid_curve.curve(),
        point: ThermalMaterialPoint {
            material: &fluid,
            context: &fluid_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    });
    let base = f.thermal();
    let mut events = base
        .cells()
        .values()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell.clone()),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: COOLANT,
            segment_id: 7,
            body_id: 7,
            material_id: MaterialId(1),
            mass_kg: 2.0,
            specific_h_j_kg: 0.5,
            geometry: CellGeometry {
                start_m: [20.0, 0.0, 0.0],
                end_m: [22.0, 0.0, 0.0],
                volume_m3: 2.0,
                surface_area_m2: 10.0,
            },
            boundary: common::boundary(),
        }),
    });
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 77,
            a: IDS[0],
            b: COOLANT,
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), f.caloric.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid_curve.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 8,
            max_links: 1,
            max_events: 9,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let bulk = MaterialNetwork {
        cells: &cells,
        links: &[],
        resistors: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 19,
                kernel_id: 74,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, &BTreeMap::from([(77, 0.0)]))
            .unwrap();
        let cooling = FixedCooling {
            prepared: &prepared,
            conductance_w_k: 0.25,
        };
        let resistors = drive_resistors(&f);
        let devices = f.devices();
        let network = CommutatedBatteryMotorHeatNetwork {
            assembly: &cooling,
            resistors: &resistors,
            devices: &devices,
            bearing_cell: IDS[0],
        };
        let steady = steady_drive(false);
        let drive = network
            .bind(&cx, &thermal, &steady.native, steady.state)
            .unwrap();
        let dt = 1.0 / 64.0;
        let next = network
            .step_candidate(
                &cx,
                &thermal,
                dt,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: dt,
                }],
                &[steady.command],
                4,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(8),
                &[0.0; 8],
            )
            .unwrap();
        let q = per_body(steady.heat_w.map(|p| dt * p), 0.0);
        // Independent two-inventory endpoint elimination:
        // Ch*dTh=Q-G*h*(Th-Tc), Cc*dTc=G*h*(Th-Tc), Ch=1,Cc=2.
        let exchange = dt * 0.25 * (0.5 + q[0]) / (1.0 + dt * 0.25 * (1.0 + 0.5));
        near(
            next.thermal().temperature_k(IDS[0]).unwrap(),
            300.0 + q[0] - exchange,
        );
        near(
            next.thermal().temperature_k(COOLANT).unwrap(),
            299.5 + exchange / 2.0,
        );
        for i in 1..7 {
            near(next.thermal().temperature_k(IDS[i]).unwrap(), 300.0 + q[i]);
        }
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            q.iter().sum(),
        );
        near(next.report().deposited_joule_j, q.iter().sum());
        near(next.report().non_electrical_heat_j, 0.0);
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.network.material.cells.len(), 8);
        near(
            next.drive().state()[0],
            steady.state[0] - dt * steady.cell_current_a,
        );
        assert!(next.thermal().temperature_k(COOLANT).unwrap() > 299.5);
    });
}

// The accepted-state check compares every source a step resolves with the bound one. The bound ones
// come from the binding's selection (SingleClaimOnly), the resolved ones from a refresh pinned to the
// claim already selected; the same sources at the same state must compare equal.
#[test]
fn a_drive_bound_under_single_claim_only_takes_its_first_step_and_another_claim_is_another_source() {
    let f = fixture(1.0, 0.0, 350.0, 100.0);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f);
    let devices = f.devices();
    let network = CommutatedBatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
        bearing_cell: IDS[0],
    };
    let steady = steady_drive(false);
    let drive = network
        .bind(&thermal, &steady.native, steady.state, || false)
        .unwrap();
    let dt = 1.0 / 64.0;
    let intervals = [CommutationInterval {
        gates: steady.gates,
        duration_s: dt,
    }];
    match network.step_candidate(
        &thermal,
        dt,
        &drive,
        &intervals,
        &[steady.command],
        4,
        None,
        power_controls(),
        cycle_controls(),
        drive_controls(7),
        &[0.0; 7],
        || false,
    ) {
        Ok(_) => {}
        Err(error) => panic!("a step from the bound state is refused: {error:?}"),
    }
    // The comparison is not vacuous: another claim value is another source.
    // (Sources refreshed by the step sit at the stepped state, so their query point differs from the
    // bound one by design; the accepted step above is the positive.)
    let other = fixture(1.0, 0.0, 350.0, 50.0);
    let other_thermal = other.thermal();
    let other_cells = other.cells(&context);
    let other_resistors = drive_resistors(&other);
    let other_devices = other.devices();
    let other_network = CommutatedBatteryMotorThermalNetwork {
        cells: &other_cells,
        links: &[],
        resistors: &other_resistors,
        devices: &other_devices,
        bearing_cell: IDS[0],
    };
    let bound_other = other_network
        .bind(&other_thermal, &steady.native, steady.state, || false)
        .unwrap();
    assert_ne!(bound_other.devices(), drive.devices());
}

// Planted negative: the drive stays bound at 300 K, the accepted heat state is another one at the same
// time, and the diode threshold depends on temperature, so the resolved sources differ in their value
// bits. The step must still refuse with the retained-source message.
#[test]
fn a_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state_at_the_same_time() {
    let f = fixture(1.0, 0.02, 350.0, 100.0);
    let thermal = f.thermal();
    let hot = f.thermal_at(320.0);
    assert!((hot.temperature_k(IDS[1]).unwrap() - 320.0).abs() < 1e-5);
    assert_eq!(hot.time_s().to_bits(), thermal.time_s().to_bits());
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f);
    let devices = f.devices();
    let network = CommutatedBatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
        bearing_cell: IDS[0],
    };
    let steady = steady_drive(false);
    let drive = network
        .bind(&thermal, &steady.native, steady.state, || false)
        .unwrap();
    let dt = 1.0 / 64.0;
    let intervals = [CommutationInterval {
        gates: steady.gates,
        duration_s: dt,
    }];
    let attempt = |at: &DepositionSimulation<'_>| {
        network
            .step_candidate(
                at,
                dt,
                &drive,
                &intervals,
                &[steady.command],
                4,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
                || false,
            )
            .map(|_| ())
            .map_err(|error| format!("{error:?}"))
    };
    assert_eq!(attempt(&thermal), Ok(()));
    let Err(text) = attempt(&hot) else {
        panic!("a step from another accepted state was taken");
    };
    assert!(
        text.contains(
            "commutated drive must retain every accepted resistance, conduction and capacitance source"
        ),
        "{text}"
    );
}
