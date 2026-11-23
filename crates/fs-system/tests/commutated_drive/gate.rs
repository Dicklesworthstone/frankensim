//! The drive's gate lane: fifteen coordinates, eight rails and four driver
//! heat ports, held to the owner built by hand, to the closed form of a gate
//! charging through its driver, and to the kernel's ledger.

use std::sync::Mutex;

use super::*;
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::PreparedGateBatteryMotorCycle;
use fs_phs::circuit::converter::commutation::capacitance::gate::GateCommutationParameters;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelPoint, GateChannelResponse, GateConductanceCurve, GateDriver,
};
use fs_ports::component::Carrier;
use fs_system::commutated_drive::{
    ChannelLaw, DRIVER_HEAT_PORTS, GateDriverResistances, GateStorage, RAIL_PORTS,
};

const CGS: [f64; 4] = [0.003, 0.004, 0.005, 0.006];
const CGD: [f64; 4] = [0.001, 0.002, 0.0015, 0.0025];
const DRIVERS: [GateDriverResistances; 4] = [
    driver(2.0, 3.0),
    driver(2.5, 2.0),
    driver(1.5, 4.0),
    driver(3.0, 1.5),
];
/// Each gate's on rail, then its off rail \[V\].
const RAILS: [f64; 8] = [3.0, -1.5, 2.8, -1.2, 3.2, -1.0, 2.6, -1.4];

const fn driver(on_resistance_ohm: f64, off_resistance_ohm: f64) -> GateDriverResistances {
    GateDriverResistances {
        on_resistance_ohm,
        off_resistance_ohm,
        maximum_current_a: 100.0,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    let bound = relative * expected.abs().max(1e-300);
    assert!(
        (actual - expected).abs() <= bound,
        "{actual:e} against {expected:e}"
    );
}

/// fs-phs's fixture channel, 0.025 S below 1 V gate-source rising to 5 S at
/// 2 V, scaled by 1 - 0.003 (T - 300 K).
fn conductance(kelvin: f64) -> GateChannelLaw {
    let scale = 1.0 - 0.003 * (kelvin - 300.0);
    let knots = [(-4.0, 0.025), (1.0, 0.025), (2.0, 5.0), (4.0, 5.0)];
    let curve = GateConductanceCurve::try_new(knots.map(|(v, g)| (v, g * scale)).to_vec());
    GateChannelLaw::new(curve.unwrap()).unwrap()
}

/// The channel of [`conductance`], refused above `upper_k`. It records each
/// temperature it is asked at.
#[derive(Debug)]
struct Channel {
    upper_k: f64,
    asked: Mutex<Vec<f64>>,
}

impl ChannelLaw for Channel {
    fn channel(&self, kelvin: f64) -> Result<GateChannelLaw, String> {
        self.asked.lock().unwrap().push(kelvin);
        if kelvin > self.upper_k {
            return Err(format!("channel: {kelvin} K is above {} K", self.upper_k));
        }
        Ok(conductance(kelvin))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.upper_k.to_bits());
    }
}

/// A channel that conducts nothing at any voltage or temperature.
#[derive(Debug, Clone, Copy)]
struct Open;

impl GateChannelResponse for Open {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-10.0, 10.0]
    }

    fn evaluate(&self, _: f64, _: f64) -> Result<GateChannelPoint, PhsError> {
        Ok(GateChannelPoint {
            current_a: 0.0,
            drain_derivative_s: 0.0,
            gate_derivative_s: 0.0,
        })
    }
}

impl ChannelLaw for Open {
    fn channel(&self, _: f64) -> Result<GateChannelLaw, String> {
        GateChannelLaw::new(*self).map_err(|error| format!("{error:?}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(0);
    }
}

/// Four channel laws of [`conductance`], each answering up to 400 K.
fn channels() -> [Arc<Channel>; 4] {
    core::array::from_fn(|_| {
        Arc::new(Channel {
            upper_k: 400.0,
            asked: Mutex::default(),
        })
    })
}

fn gates(cgd: [f64; 4], channels: [Arc<dyn ChannelLaw>; 4]) -> GateStorage {
    GateStorage {
        intrinsic_drain_charge: charge_curves(),
        gate_source_capacitance_f: CGS,
        gate_drain_capacitance_f: cgd,
        maximum_cells_per_chord: 4096,
        channels,
        drivers: DRIVERS,
    }
}

/// The owner's gate parameters on `converter`, with `channels` and the
/// drivers at `rails`.
fn gate_parameters(
    converter: CommutationParameters,
    cgd: [f64; 4],
    channels: [GateChannelLaw; 4],
    rails: &[f64],
) -> GateCommutationParameters<GateChannelLaw> {
    GateCommutationParameters {
        converter,
        intrinsic_drain_charge_curves: charge_curves(),
        gate_source_capacitance_f: CGS,
        gate_drain_capacitance_f: cgd,
        channel_conductance: channels,
        gate_drivers: core::array::from_fn(|k| GateDriver {
            on_voltage_v: rails[2 * k],
            off_voltage_v: rails[2 * k + 1],
            on_resistance_ohm: DRIVERS[k].on_resistance_ohm,
            off_resistance_ohm: DRIVERS[k].off_resistance_ohm,
            maximum_current_a: DRIVERS[k].maximum_current_a,
        }),
        maximum_cells_per_chord: 4096,
    }
}

/// The fixture's cell, motor, current and bus and node voltages, with the
/// gates at `vgs`.
fn gate_start(cgd: [f64; 4], vgs: [f64; 4]) -> [f64; 15] {
    let eleven = start(&DeviceStorage::Capacitance(COSS));
    let channels = [300.0; 4].map(conductance);
    let p = gate_parameters(parameters(100.0).converter, cgd, channels, &RAILS);
    let mut state = [0.0; 15];
    state[..3].copy_from_slice(&eleven[..3]);
    let converter = p.state(1.0, [12.0, 7.0], [-0.6, 7.6], vgs).unwrap();
    state[3..12].copy_from_slice(&converter);
    state[12..].copy_from_slice(&eleven[8..]);
    state
}

/// The drive's torques and temperatures with the rails after the torques and
/// four more temperatures, one per heat port.
fn gate_inputs(spread_k: f64) -> Vec<f64> {
    let mut inputs = inputs(spread_k);
    inputs.splice(2..2, RAILS);
    inputs.extend((11..15).map(|k| 300.0 + spread_k * f64::from(k)));
    inputs
}

#[test]
fn the_gate_lane_is_the_owners_fifteen_state_cycle_at_its_rails_and_laws() {
    let laws = channels();
    let erased = laws.clone().map(|law| law as Arc<dyn ChannelLaw>);
    let component = following(drive(gates(CGD, erased), 100.0), 400.0);
    let ports = component.ports();
    let names: Vec<_> = ports.iter().map(|port| port.id().as_str()).collect();
    let heat_ports = [HEAT_PORTS.as_slice(), &DRIVER_HEAT_PORTS, &[LOAD_PORT]].concat();
    assert_eq!(names[2..10], RAIL_PORTS);
    assert_eq!(names[10..], heat_ports);
    assert!(
        ports[2..10]
            .iter()
            .all(|port| port.carrier() == Carrier::Power(PortKind::ElectricalVoltageCurrent))
    );
    let charter = CommutatedBatteryMotorDrive::<GateStorage>::STATE_CHARTER.unwrap();
    assert!(charter.schema_grammar.starts_with("[f64; 15]"));
    let old = gate_start(CGD, [-0.5, 2.6, 2.9, -0.8]);
    let inputs = gate_inputs(7.0);
    let (torques, kelvin) = ([inputs[0], inputs[1]], &inputs[10..]);
    // The owner built by hand at the laws' answers and the rails' voltages.
    let (original, devices, commands) = owner(kelvin, torques);
    let converter = CommutationParameters {
        winding_resistance_ohm: Linear(0.1).ohms(kelvin[0]).unwrap(),
        devices,
        ..parameters(100.0).converter
    };
    let at = core::array::from_fn(|i| conductance(kelvin[1 + i]));
    let (electrical, cycle) = controls();
    let g = program();
    let mut next = [0.0; 15];
    let step = PreparedGateBatteryMotorCycle::new_gate_controlled(
        original,
        gate_parameters(converter, CGD, at, &inputs[2..10]),
        &g.intervals,
        g.substeps_per_interval,
        electrical,
        cycle,
        || false,
    )
    .unwrap()
    .step_interval_commands_controlled(&old, &commands, &mut next, || false)
    .unwrap();
    let account = step.converter.gate.as_ref().unwrap();
    // The switched gates draw on both rails; the others on their off rails.
    let q = account.rail_charge_c;
    assert!(q[1][0] > 0.0 && q[2][0] > 0.0 && q.iter().all(|q| q[1] < 0.0));
    assert!(q[0][0] == 0.0 && q[3][0] == 0.0, "{q:?}");
    let mut heats = port_heats(step.heat_joules(), step.load_work_j).to_vec();
    heats.splice(10..10, account.resistor_heat_j);
    let expected: Vec<u64> = step
        .mean_shaft_speed_rad_s
        .into_iter()
        .chain(q.as_flattened().iter().map(|q| q / PERIOD))
        .chain(heats.iter().map(|q| -q / PERIOD))
        .map(f64::to_bits)
        .collect();
    with_cx(|cx| {
        let (mut new, mut outputs) = ([0.0; 15], [0.0; 25]);
        let proposal = component
            .propose(&old, PERIOD, &inputs, &mut new, &mut outputs, cx)
            .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the drive shortened the step");
        };
        assert_eq!(new.map(f64::to_bits), next.map(f64::to_bits));
        assert_eq!(outputs.map(f64::to_bits).to_vec(), expected);
        // Each channel law was asked once, at its own switch's heat port.
        for (i, law) in laws.iter().enumerate() {
            assert_eq!(*law.asked.lock().unwrap(), [kelvin[1 + i]]);
        }
        let stored =
            component.stored_energy_j(&new).unwrap() - component.stored_energy_j(&old).unwrap();
        assert!((receipt.stored_change_j - stored).abs() <= receipt.energy_tolerance_j);
        close(
            receipt.stored_change_j,
            step.total_storage_change_j(),
            1e-12,
        );
        let entropy: f64 = heats.iter().zip(kelvin).map(|(q, t)| q / t).sum();
        close(receipt.entropy_generated_j_per_k, entropy, 1e-15);
        let case = Case {
            state: &old,
            stale: &old.map(|q| 1.01 * q),
            dt_s: PERIOD,
            inputs: &inputs,
        };
        check_component(&component, &case, cx).expect("component conformance");
    });
}

#[test]
fn a_refused_gate_step_leaves_nothing_behind_and_the_identity_names_the_lane() {
    let laws = channels().map(|law| law as Arc<dyn ChannelLaw>);
    // The switch laws answer to 450 K, so the channel's own refusal shows.
    let component = following(drive(gates(CGD, laws.clone()), 100.0), 450.0);
    let old = gate_start(CGD, [-0.5, 2.6, 2.9, -0.8]);
    let inputs = gate_inputs(7.0);
    let refuse = |inputs: &[f64]| refused(&component, &old, PERIOD, inputs);
    // A non-finite rail; a zero driver temperature.
    for (index, value) in [(3, f64::NAN), (21, 0.0)] {
        let mut bad = inputs.clone();
        bad[index] = value;
        let error = refuse(&bad);
        assert!(matches!(error, CommutatedDriveError::Invalid { .. }));
    }
    // An on rail driving the input low gate past its driver's current limit
    // is the owner's refusal.
    let mut high = inputs.clone();
    high[4] = 1000.0;
    let error = refuse(&high);
    let text = error.to_string();
    assert!(matches!(error, CommutatedDriveError::Step(_)) && text.contains("gate-driver current"));
    // The input low channel above its law's range.
    let mut hot = inputs.clone();
    hot[12] = 410.0;
    let error = refuse(&hot);
    let CommutatedDriveError::Law { port, reason } = error else {
        panic!("{error}");
    };
    assert!(port == "heat_input_low" && reason.starts_with("channel"));
    // The codec keeps fifteen words; the identity names the lane.
    let mut words = Vec::new();
    assert!(component.encode_state(&old, &mut words));
    assert_eq!(component.decode_state(&words), Some(old));
    assert_eq!(component.decode_state(&words[..11]), None);
    let identity = |storage| {
        let mut words = Vec::new();
        drive(storage, 100.0).identity_words(&mut words);
        words
    };
    let open = || [(); 4].map(|()| Arc::new(Open) as Arc<dyn ChannelLaw>);
    let gated = identity(gates(CGD, open()));
    let mut other = Vec::new();
    drive(DeviceStorage::Capacitance(COSS), 100.0).identity_words(&mut other);
    let part = gated.iter().zip(&other).position(|(a, b)| a != b).unwrap();
    assert_eq!((gated[part], other[part]), (3, 0));
    assert_ne!(gated, identity(gates([0.0; 4], open())));
    let mut slower = gates(CGD, open());
    slower.drivers[2].off_resistance_ohm = 5.0;
    assert_ne!(gated, identity(slower));
    assert_ne!(gated, identity(gates(CGD, laws)));
}

/// Every gate commanded `on` for half a period, then off; no load.
fn edge_program(on: [bool; 4]) -> GateProgram {
    GateProgram {
        intervals: [on, [false; 4]]
            .map(|on| CommutationInterval {
                gates: DeviceGateState { on },
                duration_s: 0.5 * PERIOD,
            })
            .to_vec(),
        auxiliary_loads: vec![None; 2],
        substeps_per_interval: 4,
    }
}

#[test]
fn a_gate_without_miller_coupling_charges_in_closed_form_and_moves_nothing_else() {
    let (electrical, cycle) = controls();
    let open = || [(); 4].map(|()| Arc::new(Open) as Arc<dyn ChannelLaw>);
    let built = |on| {
        let storage = gates([0.0; 4], open());
        let drive = CommutatedBatteryMotorDrive::new(
            parameters(100.0),
            storage,
            edge_program(on),
            electrical,
            cycle,
            1e-9,
        );
        following(drive.unwrap(), 400.0)
    };
    let old = gate_start([0.0; 4], [0.0; 4]);
    let all = gate_inputs(5.0);
    let inputs = &all[..24];
    let run = |on| {
        with_cx(|cx| {
            let (mut new, mut outputs) = ([0.0; 15], [0.0; 24]);
            built(on)
                .propose(&old, PERIOD, inputs, &mut new, &mut outputs, cx)
                .unwrap();
            (new, outputs)
        })
    };
    let (next, outputs) = run([true; 4]);
    let h = 0.5 * PERIOD / 4.0;
    let mut worst = 0.0_f64;
    for k in 0..4 {
        let c = CGS[k];
        let [on, off] = [RAILS[2 * k], RAILS[2 * k + 1]];
        // A constant rail through one resistance into Cgs: the owner's
        // implicit midpoint takes v - V by (2RC - h)/(2RC + h) per substep.
        let edge = |v: f64, rail: f64, ohm: f64| {
            rail + (v - rail) * ((2.0 * ohm * c - h) / (2.0 * ohm * c + h)).powi(4)
        };
        let v1 = edge(0.0, on, DRIVERS[k].on_resistance_ohm);
        let v2 = edge(v1, off, DRIVERS[k].off_resistance_ohm);
        let charge = [c * v1, c * (v2 - v1)];
        // Each rail's work less the gate's stored change.
        let heat = on * charge[0] + off * charge[1] - 0.5 * c * v2 * v2;
        let measured = [
            next[8 + k] / c,
            outputs[2 + 2 * k] * PERIOD,
            outputs[3 + 2 * k] * PERIOD,
            -outputs[20 + k] * PERIOD,
        ];
        for (value, exact) in measured.into_iter().zip([v2, charge[0], charge[1], heat]) {
            worst = worst.max((value - exact).abs() / exact.abs());
            close(value, exact, 1e-12);
        }
    }
    // With the channels open and no Miller coupling, the gates reach nothing
    // else: every other coordinate, both shafts and the ten drive heat ports
    // are those of a program whose gates never switch.
    let (still, still_outputs) = run([false; 4]);
    for i in (0..8).chain(12..15) {
        assert_eq!(next[i].to_bits(), still[i].to_bits(), "coordinate {i}");
    }
    for i in (0..2).chain(10..20) {
        assert_eq!(
            outputs[i].to_bits(),
            still_outputs[i].to_bits(),
            "output {i}"
        );
    }
    assert!((8..12).all(|i| next[i].to_bits() != still[i].to_bits()));
    println!("gate RC closed form: {worst:.1e} relative");
}

#[test]
fn the_kernels_ledger_closes_with_a_source_on_every_rail() {
    let mut builder = SystemBuilder::new();
    let laws = channels().map(|law| law as Arc<dyn ChannelLaw>);
    let component = following(drive(gates(CGD, laws), 100.0), 400.0);
    let old = gate_start(CGD, [-0.5, 2.6, 2.9, -0.8]);
    let node = builder.add("drive", component, old).unwrap();
    let efforts = [("rotor", -0.05), ("stator", 0.03)]
        .map(|(name, value)| (name, value, PortKind::RotationalTorqueAngularVelocity));
    let rails = RAIL_PORTS
        .into_iter()
        .zip(RAILS)
        .map(|(name, value)| (name, value, PortKind::ElectricalVoltageCurrent));
    for (name, value, kind) in efforts.into_iter().chain(rails) {
        let source = Force::in_domain(kind, Waveform::constant(value).unwrap());
        let source = builder.add(name, source, SourceAccount::default()).unwrap();
        let (a, b) = (
            builder.port(source, "p").unwrap(),
            builder.port(node, name).unwrap(),
        );
        builder.bond(a, b).unwrap();
    }
    let names = [HEAT_PORTS.as_slice(), &DRIVER_HEAT_PORTS, &[LOAD_PORT]].concat();
    for name in &names {
        let body = builder
            .add(
                &format!("body {name}"),
                HeatCapacity::new(1e-2).unwrap(),
                300.0,
            )
            .unwrap();
        let (a, b) = (
            builder.port(body, "heat").unwrap(),
            builder.port(node, name).unwrap(),
        );
        builder.common_effort(a, &[b]).unwrap();
    }
    let mut system = builder.seal().unwrap();
    for name in &names {
        let port = system.port(node, name).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    with_cx(|cx| {
        for _ in 0..3 {
            system.step(PERIOD, cx).unwrap();
        }
    });
    let after = *system.state::<[f64; 15]>(node).unwrap();
    assert_ne!(after.map(f64::to_bits), old.map(f64::to_bits));
    let residual = system.ledger().energy_residual_j();
    assert!(residual.abs() <= 1e-15, "ledger residual {residual}");
    println!("gate lane on the kernel: ledger {residual:.1e} J");
}
