//! The battery-fed commutated motor drive as one Component: a step is its
//! owner's gate program at the laws' answers, the kernel steps it at the
//! temperatures its heat ports receive, and every refusal leaves the outputs and the candidate
//! untouched.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::SolveBudget;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    PreparedCapacitiveBatteryMotorCycle, PreparedChargeBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::charge::ChargeCommutationParameters;
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState, DeviceOperatingBounds,
    SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, check_component};
use fs_ports::component::{Component, Proposal};
use fs_system::SystemBuilder;
use fs_system::circuit::ResistanceLaw;
use fs_system::commutated_drive::{
    CommutatedBatteryMotorDrive, CommutatedDriveError, DeviceStorage, DriveParameters,
    DriveResistor, GateProgram, HEAT_PORTS, LOAD_PORT, SwitchLaw, SwitchPosition,
};
use fs_system::elements::{Force, HeatCapacity, SourceAccount, Waveform};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 5,
            kernel_id: 61,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

const COSS: [f64; 4] = [0.025, 0.04, 0.035, 0.05];
const OCV: [(f64, f64); 3] = [(0.0, 11.0), (0.5, 12.0), (1.0, 12.6)];
const PERIOD: f64 = 1.0 / 256.0;

fn switch(forward_voltage_v: f64, channel_limit_a: f64) -> SwitchDevice {
    SwitchDevice {
        on_resistance_ohm: 0.1,
        off_resistance_ohm: 100.0,
        body_diode: BodyDiode {
            forward_voltage_v,
            differential_resistance_ohm: 0.2,
        },
        bounds: DeviceOperatingBounds {
            maximum_reverse_voltage_v: 100.0,
            maximum_channel_current_a: channel_limit_a,
            maximum_diode_current_a: 100.0,
        },
    }
}

fn parameters(channel_limit_a: f64) -> DriveParameters {
    DriveParameters {
        capacity_c: 100.0,
        soc_voltage: OCV.to_vec(),
        battery: BatteryParameters {
            series_resistance_ohm: 0.1,
            polarization_resistance_ohm: 0.1,
            polarization_capacitance_f: 2.0,
            bus_capacitance_f: 0.5,
        },
        lead_resistance_ohm: 0.1,
        converter: CommutationParameters {
            inductance_h: 0.3,
            bus_capacitance_f: [0.2, 0.4],
            winding_resistance_ohm: 0.1,
            devices: [switch(0.5, channel_limit_a); 4],
        },
        motor: DcMotorParameters {
            inductance_h: 0.05,
            winding_resistance_ohm: 0.1,
            rotor_inertia_kg_m2: 0.4,
            stator_inertia_kg_m2: 0.8,
            coupling_nm_per_a: 0.2,
            bearing_drag_nm_s_per_rad: 0.0003,
        },
    }
}

/// A dead time, then input low and output high on, with an auxiliary load
/// throughout. The fixture's nodes start where the input low and output
/// high diodes conduct, so both carry current in the dead time.
fn program() -> GateProgram {
    let load = RegulatedPowerLoad::try_new(0.5, 10.0, 0.25, 100.0).unwrap();
    let on = [([false; 4], 0.1), ([false, true, true, false], 0.9)];
    GateProgram {
        intervals: on
            .map(|(on, fraction)| CommutationInterval {
                gates: DeviceGateState { on },
                duration_s: fraction * PERIOD,
            })
            .to_vec(),
        auxiliary_loads: vec![Some(load); 2],
        substeps_per_interval: 3,
    }
}

fn controls() -> (PowerLoadControls, SwitchCycleControls) {
    (
        PowerLoadControls {
            solve: SolveBudget {
                relative: 1e-13,
                absolute: 1e-12,
                max_iterations: 40,
            },
            backtracking_halvings: 10,
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-12,
        },
        SwitchCycleControls {
            max_intervals: 64,
            absolute_energy_tolerance_j: 1e-9,
            relative_energy_tolerance: 1e-11,
        },
    )
}

fn charge_curves() -> [CapacitorChargeCurve; 4] {
    let volts = [-100.0, -5.0, -0.5, 0.0, 0.5, 2.0, 20.0, 100.0];
    let charge = [-50.0, -5.0, -1.0, 0.0, 1.0, 2.5, 7.0, 15.0];
    COSS.map(|scale| {
        CapacitorChargeCurve::try_new(charge.iter().map(|q| scale * q).zip(volts).collect())
            .unwrap()
    })
}

fn storages() -> [DeviceStorage; 2] {
    [
        DeviceStorage::Capacitance(COSS),
        DeviceStorage::Charge {
            curves: charge_curves(),
            maximum_cells_per_chord: 4096,
        },
    ]
}

fn drive(storage: DeviceStorage, channel_limit_a: f64) -> CommutatedBatteryMotorDrive {
    let (electrical, cycle) = controls();
    CommutatedBatteryMotorDrive::new(
        parameters(channel_limit_a),
        storage,
        program(),
        electrical,
        cycle,
        1e-9,
    )
    .unwrap()
}

/// A resistance that rises linearly with temperature.
#[derive(Debug)]
struct Linear(f64);

impl ResistanceLaw for Linear {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        Ok(self.0 * (1.0 + 0.004 * (kelvin - 300.0)))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.0.to_bits());
    }
}

/// The fixture's switch with a body diode whose threshold falls 2 mV/K,
/// refused above `upper_k`.
#[derive(Debug)]
struct Threshold {
    upper_k: f64,
}

impl SwitchLaw for Threshold {
    fn switch(&self, kelvin: f64) -> Result<SwitchDevice, String> {
        if kelvin > self.upper_k {
            return Err(format!("{kelvin} K is above {} K", self.upper_k));
        }
        Ok(switch(0.5 - 0.002 * (kelvin - 300.0), 100.0))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.upper_k.to_bits());
    }
}

const RESISTORS: [DriveResistor; 5] = [
    DriveResistor::ConverterWinding,
    DriveResistor::BatterySeries,
    DriveResistor::BatteryPolarization,
    DriveResistor::Lead,
    DriveResistor::MotorWinding,
];
const SWITCHES: [SwitchPosition; 4] = [
    SwitchPosition::InputHigh,
    SwitchPosition::InputLow,
    SwitchPosition::OutputHigh,
    SwitchPosition::OutputLow,
];

fn following(mut drive: CommutatedBatteryMotorDrive, upper_k: f64) -> CommutatedBatteryMotorDrive {
    for resistor in RESISTORS {
        drive = drive.following(resistor, Arc::new(Linear(0.1))).unwrap();
    }
    for position in SWITCHES {
        drive = drive
            .switch_following(position, Arc::new(Threshold { upper_k }))
            .unwrap();
    }
    drive
}

fn start(storage: &DeviceStorage) -> [f64; 11] {
    let p = parameters(100.0);
    let curve = BatteryOcvCurve::new(p.capacity_c, &OCV).unwrap();
    let battery = BatterySource::new(curve, p.battery).unwrap();
    let cell = battery.state(0.8, 0.05, 12.2).unwrap();
    let (buses, nodes) = ([12.0, 7.0], [-0.6, 7.6]);
    let converter = match storage {
        DeviceStorage::Capacitance(c) => CapacitiveCommutationParameters {
            converter: p.converter,
            output_capacitance_f: *c,
        }
        .state(1.0, buses, nodes),
        DeviceStorage::Charge {
            curves,
            maximum_cells_per_chord,
        } => ChargeCommutationParameters {
            converter: p.converter,
            output_charge_curves: curves.clone(),
            maximum_cells_per_chord: *maximum_cells_per_chord,
        }
        .state(1.0, buses, nodes),
    }
    .unwrap();
    let motor = DcMotor::new(p.motor)
        .unwrap()
        .state(0.8, [30.0, 2.0])
        .unwrap();
    let mut state = [0.0; 11];
    state[..3].copy_from_slice(&cell);
    state[3..8].copy_from_slice(&converter);
    state[8..].copy_from_slice(&motor);
    state
}

/// Torques on the rotor and the stator, then one temperature per heat port.
fn inputs(spread_k: f64) -> Vec<f64> {
    let mut inputs = vec![-0.05, 0.03];
    inputs.extend((0..11).map(|k| 300.0 + spread_k * f64::from(k)));
    inputs
}

/// The owner's own step at what the laws answer at `kelvin`, built here
/// without the Component.
fn native(
    storage: &DeviceStorage,
    old: &[f64; 11],
    torques: [f64; 2],
    kelvin: &[f64],
) -> ([f64; 11], [f64; 2], [f64; 11], [f64; 14]) {
    let ohms = |k: usize| Linear(0.1).ohms(kelvin[k]).unwrap();
    let mut p = parameters(100.0);
    p.converter.winding_resistance_ohm = ohms(0);
    for i in 0..4 {
        p.converter.devices[i] = Threshold { upper_k: 400.0 }.switch(kelvin[1 + i]).unwrap();
    }
    let on = p.converter.devices.map(|d| d.on_resistance_ohm);
    p.battery.series_resistance_ohm = ohms(5);
    p.battery.polarization_resistance_ohm = ohms(6);
    p.motor.winding_resistance_ohm = ohms(8);
    let original: BatteryMotorDrive = AveragedBuckBoost::new(
        BuckBoostParameters {
            inductance_h: 0.3,
            bus_capacitance_f: [0.2, 0.4],
            winding_resistance_ohm: ohms(0),
            input_on_resistance_ohm: [on[0], on[1]],
            output_on_resistance_ohm: [on[2], on[3]],
        },
        BridgeModulation {
            input_high_fraction: 0.0,
            output_high_fraction: 0.0,
        },
    )
    .unwrap()
    .with_battery_source(
        BatterySource::new(BatteryOcvCurve::new(100.0, &OCV).unwrap(), p.battery).unwrap(),
        ohms(7),
    )
    .unwrap()
    .with_dc_motor(DcMotor::new(p.motor).unwrap())
    .unwrap();
    let (electrical, cycle) = controls();
    let g = program();
    let commands: Vec<_> = g
        .auxiliary_loads
        .iter()
        .map(|&auxiliary_load| BatteryMotorSwitchCommand {
            shaft_torque_nm: torques,
            auxiliary_load,
        })
        .collect();
    let mut next = [0.0; 11];
    let devices = p.converter.devices;
    let (n, k) = (&g.intervals, g.substeps_per_interval);
    let step = match storage {
        DeviceStorage::Capacitance(c) => PreparedCapacitiveBatteryMotorCycle::new_controlled(
            original,
            devices,
            *c,
            n,
            k,
            electrical,
            cycle,
            || false,
        )
        .unwrap()
        .step_interval_commands_controlled(old, &commands, &mut next, || false),
        DeviceStorage::Charge {
            curves,
            maximum_cells_per_chord,
        } => PreparedChargeBatteryMotorCycle::new_charge_controlled(
            original,
            devices,
            curves.clone(),
            *maximum_cells_per_chord,
            n,
            k,
            electrical,
            cycle,
            || false,
        )
        .unwrap()
        .step_interval_commands_controlled(old, &commands, &mut next, || false),
    }
    .unwrap();
    let h = step.heat_joules();
    let heats = [
        h[0],
        h[1] + h[10],
        h[2] + h[11],
        h[3] + h[12],
        h[4] + h[13],
        h[5],
        h[6],
        h[7],
        h[8],
        h[9],
        step.load_work_j,
    ];
    (next, step.mean_shaft_speed_rad_s, heats, h)
}

#[test]
fn a_step_is_the_owners_gate_program_at_what_the_laws_answer() {
    for storage in storages() {
        let old = start(&storage);
        let component = following(drive(storage.clone(), 100.0), 400.0);
        assert_eq!(component.period_s().to_bits(), PERIOD.to_bits());
        assert_eq!(component.ports().len(), 13);
        assert_eq!(component.ports()[12].id().as_str(), LOAD_PORT);
        let inputs = inputs(7.0);
        let (next, speeds, heats, owner) =
            native(&storage, &old, [inputs[0], inputs[1]], &inputs[2..]);
        // Both diodes the program puts in conduction carry heat of their own.
        assert!(
            owner[11] > 0.0 && owner[12] > 0.0 && owner[11].to_bits() != owner[12].to_bits(),
            "{owner:?}"
        );
        assert!(
            heats[1..5].iter().all(|q| *q > 0.0) && heats[10] > 0.0,
            "{heats:?}"
        );
        with_cx(|cx| {
            let mut new = [0.0; 11];
            let mut outputs = [0.0; 13];
            let proposal = component
                .propose(&old, PERIOD, &inputs, &mut new, &mut outputs, cx)
                .unwrap();
            let Proposal::Step(receipt) = proposal else {
                panic!("the drive shortened the step");
            };
            assert_eq!(new.map(f64::to_bits), next.map(f64::to_bits));
            assert_eq!(outputs[0].to_bits(), speeds[0].to_bits());
            assert_eq!(outputs[1].to_bits(), speeds[1].to_bits());
            for (output, heat) in outputs[2..].iter().zip(heats) {
                assert_eq!(output.to_bits(), (-heat / PERIOD).to_bits());
            }
            let stored = component.stored_energy_j(&next).unwrap()
                - component.stored_energy_j(&old).unwrap();
            assert!((receipt.stored_change_j - stored).abs() <= receipt.energy_tolerance_j);
            let entropy: f64 = heats.iter().zip(&inputs[2..]).map(|(q, t)| q / t).sum();
            assert!((receipt.entropy_generated_j_per_k - entropy).abs() <= 1e-15 * entropy);
            let case = Case {
                state: &old,
                stale: &start(&DeviceStorage::Capacitance(COSS)).map(|q| 1.01 * q),
                dt_s: PERIOD,
                inputs: &inputs,
            };
            check_component(&component, &case, cx).expect("component conformance");
        });
    }
}

#[test]
fn the_kernel_steps_the_drive_at_its_bodies_mean_temperatures() {
    let capacities = [
        2e-3, 1e-3, 1e-3, 1e-3, 1e-3, 5e-3, 5e-3, 5e-3, 5e-3, 1.0, 2e-3,
    ];
    for storage in storages() {
        let old = start(&storage);
        let mut builder = SystemBuilder::new();
        let node = builder
            .add(
                "drive",
                following(drive(storage.clone(), 100.0), 400.0),
                old,
            )
            .unwrap();
        for (shaft, torque) in [("rotor", -0.05), ("stator", 0.03)] {
            let kind = PortKind::RotationalTorqueAngularVelocity;
            let source = Force::in_domain(kind, Waveform::constant(torque).unwrap());
            let source = builder
                .add(shaft, source, SourceAccount::default())
                .unwrap();
            let (a, b) = (
                builder.port(source, "p").unwrap(),
                builder.port(node, shaft).unwrap(),
            );
            builder.bond(a, b).unwrap();
        }
        let names: Vec<_> = HEAT_PORTS.into_iter().chain([LOAD_PORT]).collect();
        let mut bodies = Vec::new();
        for (name, capacity) in names.iter().zip(capacities) {
            let body = HeatCapacity::new(capacity).unwrap();
            let body = builder.add(&format!("body {name}"), body, 300.0).unwrap();
            let (a, b) = (
                builder.port(body, "heat").unwrap(),
                builder.port(node, name).unwrap(),
            );
            builder.common_effort(a, &[b]).unwrap();
            bodies.push(body);
        }
        let mut system = builder.seal().unwrap();
        for name in &names {
            let port = system.port(node, name).unwrap();
            assert!(system.guess(port, &[300.0]));
        }
        let check = following(drive(storage.clone(), 100.0), 400.0);
        let mut kelvin = [300.0; 11];
        let mut worst = [0.0_f64; 2];
        with_cx(|cx| {
            for _ in 0..6 {
                let before = *system.state::<[f64; 11]>(node).unwrap();
                system.step(PERIOD, cx).unwrap();
                let after = *system.state::<[f64; 11]>(node).unwrap();
                let end: Vec<f64> = bodies
                    .iter()
                    .map(|&body| *system.state::<f64>(body).unwrap())
                    .collect();
                // The drive proposed at the temperatures the bodies return,
                // their means over the step, takes the kernel's step and
                // gives each body what it received.
                let mut given = vec![-0.05, 0.03];
                given.extend(end.iter().zip(kelvin).map(|(b, a)| f64::midpoint(a, *b)));
                let (mut new, mut outputs) = ([0.0; 11], [0.0; 13]);
                check
                    .propose(&before, PERIOD, &given, &mut new, &mut outputs, cx)
                    .unwrap();
                for (i, (a, b)) in new.iter().zip(after).enumerate() {
                    let scale = a.abs().max(1.0);
                    worst[0] = worst[0].max((a - b).abs() / scale);
                    assert!(
                        (a - b).abs() <= 1.5e-15 * scale,
                        "coordinate {i}: {a} vs {b}"
                    );
                }
                for (k, rate) in outputs[2..].iter().enumerate() {
                    let received = capacities[k] * (end[k] - kelvin[k]);
                    let heat = -rate * PERIOD;
                    worst[1] = worst[1].max((heat - received).abs());
                    assert!(
                        (heat - received).abs() <= 2.5e-13,
                        "{}: {heat} vs {received}",
                        names[k]
                    );
                }
                kelvin.copy_from_slice(&end);
            }
        });
        assert!(
            kelvin[1] > 330.0,
            "the switches heat their bodies: {kelvin:?}"
        );
        let residual = system.ledger().energy_residual_j();
        assert!(residual.abs() <= 1e-16, "ledger residual {residual}");
        println!(
            "mean-temperature step: coordinates {:.1e}, heat {:.1e} J, ledger {residual:.1e} J, \
             input high {:.2} K",
            worst[0], worst[1], kelvin[1]
        );
    }
}

fn untouched(values: &[f64], sentinel: f64) -> bool {
    values
        .iter()
        .all(|value| value.to_bits() == sentinel.to_bits())
}

fn refused(
    drive: &CommutatedBatteryMotorDrive,
    old: &[f64; 11],
    dt_s: f64,
    inputs: &[f64],
) -> CommutatedDriveError {
    with_cx(|cx| {
        let mut new = [7.0; 11];
        let mut outputs = vec![9.0; inputs.len()];
        let error = drive
            .propose(old, dt_s, inputs, &mut new, &mut outputs, cx)
            .unwrap_err();
        assert!(untouched(&new, 7.0) && untouched(&outputs, 9.0));
        error
    })
}

#[test]
fn refusals_leave_nothing_behind() {
    let (electrical, cycle) = controls();
    let build = |program: GateProgram, tolerance: f64| {
        CommutatedBatteryMotorDrive::new(
            parameters(100.0),
            DeviceStorage::Capacitance(COSS),
            program,
            electrical,
            cycle,
            tolerance,
        )
    };
    let mut unpaired = program();
    unpaired.auxiliary_loads.pop();
    let mut empty = program();
    empty.intervals.clear();
    empty.auxiliary_loads.clear();
    assert!(matches!(
        build(program(), 0.0),
        Err(CommutatedDriveError::Invalid { .. })
    ));
    assert!(matches!(
        build(unpaired, 1e-9),
        Err(CommutatedDriveError::Invalid { .. })
    ));
    assert!(build(empty, 1e-9).is_err());
    let twice = following(drive(DeviceStorage::Capacitance(COSS), 100.0), 400.0);
    assert!(
        twice
            .following(DriveResistor::Lead, Arc::new(Linear(0.1)))
            .is_err()
    );
    let twice = following(drive(DeviceStorage::Capacitance(COSS), 100.0), 400.0);
    let law = Arc::new(Threshold { upper_k: 400.0 });
    assert!(
        twice
            .switch_following(SwitchPosition::InputLow, law)
            .is_err()
    );

    for storage in storages() {
        let old = start(&storage);
        let good = inputs(7.0);
        let component = following(drive(storage.clone(), 100.0), 350.0);
        let error = refused(&component, &old, 0.5 * PERIOD, &good);
        assert!(matches!(error, CommutatedDriveError::Invalid { .. }));
        for (index, value) in [(0, f64::NAN), (5, 0.0), (12, f64::INFINITY)] {
            let mut bad = good.clone();
            bad[index] = value;
            let error = refused(&component, &old, PERIOD, &bad);
            assert!(matches!(error, CommutatedDriveError::Invalid { .. }));
        }
        // The input low switch's package above its law's validity.
        let mut hot = good.clone();
        hot[4] = 360.0;
        let error = refused(&component, &old, PERIOD, &hot);
        assert!(
            matches!(
                error,
                CommutatedDriveError::Law {
                    port: "heat_input_low",
                    ..
                }
            ),
            "{error}"
        );
        // A current above the channels' operating limit is the owner's
        // refusal, unchanged.
        let error = refused(&drive(storage.clone(), 0.05), &old, PERIOD, &good);
        assert!(matches!(error, CommutatedDriveError::Step(_)), "{error}");
        // Cancellation before or during the step.
        let gate = CancelGate::new_clock_free();
        gate.request();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let key = StreamKey {
                seed: 5,
                kernel_id: 61,
                tile: 0,
                iteration: 0,
            };
            let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
            let (mut new, mut outputs) = ([7.0; 11], [9.0; 13]);
            let error = component
                .propose(&old, PERIOD, &good, &mut new, &mut outputs, &cx)
                .unwrap_err();
            assert!(matches!(error, CommutatedDriveError::Cancelled));
            assert!(untouched(&new, 7.0) && untouched(&outputs, 9.0));
        });
    }
}

#[test]
fn the_codec_and_the_identity_tell_states_and_drives_apart() {
    for storage in storages() {
        let old = start(&storage);
        let component = following(drive(storage, 100.0), 350.0);
        // The codec keeps a state and refuses one outside the laws.
        let mut words = Vec::new();
        assert!(component.encode_state(&old, &mut words));
        assert_eq!(component.decode_state(&words), Some(old));
        assert_eq!(component.decode_state(&words[..10]), None);
        let mut outside = old;
        outside[0] = 2.0 * 100.0;
        assert!(component.stored_energy_j(&outside).is_err());
        assert!(!component.encode_state(&outside, &mut Vec::new()));
        words[0] = outside[0].to_bits();
        assert_eq!(component.decode_state(&words), None);
    }
    // The identity holds the storage law and every law followed.
    let identity = |drive: &CommutatedBatteryMotorDrive| {
        let mut words = Vec::new();
        drive.identity_words(&mut words);
        words
    };
    let [capacitance, charge] = storages().map(|storage| drive(storage, 100.0));
    let plain = identity(&capacitance);
    assert_ne!(plain, identity(&charge));
    let lead = capacitance
        .following(DriveResistor::Lead, Arc::new(Linear(0.1)))
        .unwrap();
    assert_ne!(plain, identity(&lead));
    assert_ne!(
        identity(&lead),
        identity(&following(
            drive(DeviceStorage::Capacitance(COSS), 100.0),
            400.0
        ))
    );
}

/// The test program with only the output high switch on after the dead time.
fn other_program() -> GateProgram {
    let mut program = program();
    program.intervals[1].gates.on = [false, false, true, false];
    program
}

/// One proposal of one period; the new state and the thirteen outputs, or
/// the refusal, which leaves both untouched.
fn step(
    drive: &CommutatedBatteryMotorDrive,
    old: &[f64; 11],
    inputs: &[f64],
) -> Result<([f64; 11], [f64; 13]), CommutatedDriveError> {
    with_cx(|cx| {
        let (mut new, mut outputs) = ([7.0; 11], [9.0; 13]);
        let result = drive.propose(old, PERIOD, inputs, &mut new, &mut outputs, cx);
        if result.is_err() {
            assert!(untouched(&new, 7.0) && untouched(&outputs, 9.0));
        }
        result.map(|_| (new, outputs))
    })
}

#[test]
fn a_program_signal_chooses_each_steps_program() {
    let (electrical, cycle) = controls();
    for storage in storages() {
        let old = start(&storage);
        let built = |program| {
            let drive = CommutatedBatteryMotorDrive::new(
                parameters(100.0),
                storage.clone(),
                program,
                electrical,
                cycle,
                1e-9,
            );
            following(drive.unwrap(), 400.0)
        };
        let scheduled = built(program()).scheduled(vec![other_program()]).unwrap();
        assert_eq!(scheduled.ports().len(), 14);
        assert_eq!(scheduled.ports()[13].id().as_str(), "program");
        // Each index runs what the drive built with that program alone runs.
        let given = inputs(5.0);
        let mut ends = Vec::new();
        for (index, program) in [0.0, 1.0].into_iter().zip([program(), other_program()]) {
            let signalled = [given.as_slice(), &[index]].concat();
            let (new, outputs) = step(&scheduled, &old, &signalled).unwrap();
            let (alone, alone_outputs) = step(&built(program), &old, &given).unwrap();
            assert_eq!(new.map(f64::to_bits), alone.map(f64::to_bits));
            assert_eq!(outputs.map(f64::to_bits), alone_outputs.map(f64::to_bits));
            ends.push(new);
        }
        assert_ne!(ends[0].map(f64::to_bits), ends[1].map(f64::to_bits));
        // A signal that names no program is refused.
        for index in [0.5, 2.0, -1.0] {
            let signalled = [given.as_slice(), &[index]].concat();
            let error = step(&scheduled, &old, &signalled).unwrap_err();
            assert!(error.to_string().contains("program signal"), "{error}");
        }
        with_cx(|cx| {
            let case = Case {
                state: &old,
                stale: &start(&DeviceStorage::Capacitance(COSS)).map(|q| 1.01 * q),
                dt_s: PERIOD,
                inputs: &[given.as_slice(), &[1.0]].concat(),
            };
            check_component(&scheduled, &case, cx).expect("component conformance");
        });
        // Programs are given once, nonempty, each with an entry per interval
        // and the drive's period.
        let refuses = |programs| built(program()).scheduled(programs).is_err();
        let mut longer = other_program();
        longer.intervals[1].duration_s *= 2.0;
        let mut unloaded = other_program();
        unloaded.auxiliary_loads.pop();
        assert!(refuses(Vec::new()) && refuses(vec![longer]) && refuses(vec![unloaded]));
        assert!(scheduled.scheduled(vec![other_program()]).is_err());
        // The identity holds the further programs.
        let identity = |drive: CommutatedBatteryMotorDrive| {
            let mut words = Vec::new();
            drive.identity_words(&mut words);
            words
        };
        let plain = identity(built(program()));
        let other = identity(built(program()).scheduled(vec![other_program()]).unwrap());
        let same = identity(built(program()).scheduled(vec![program()]).unwrap());
        assert!(plain != other && plain != same && other != same);
    }
}
